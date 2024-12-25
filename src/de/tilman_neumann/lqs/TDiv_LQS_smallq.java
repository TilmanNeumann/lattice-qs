/*
 * lattics-qs (LQS) is a study of a lattice quadratic sieve proposed by R.D.Silverman in https://www.mersenneforum.org/showthread.php?t=14080.
 * Copyright (C) 2018 Tilman Neumann - tilman.neumann@web.de
 *
 * This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program;
 * if not, see <http://www.gnu.org/licenses/>.
 */
package de.tilman_neumann.lqs;

import static de.tilman_neumann.jml.factor.base.GlobalFactoringOptions.*;
import static de.tilman_neumann.jml.base.BigIntConstants.*;
import static org.junit.Assert.*;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

import de.tilman_neumann.jml.base.UnsignedBigInt;
import de.tilman_neumann.jml.factor.base.SortedIntegerArray;
import de.tilman_neumann.jml.factor.base.SortedLongArray;
import de.tilman_neumann.jml.factor.base.congruence.AQPair;
import de.tilman_neumann.jml.factor.base.congruence.AQPairFactory;
import de.tilman_neumann.jml.factor.siqs.SIQS_Small;
import de.tilman_neumann.jml.factor.siqs.poly.SIQSPolyGenerator;
import de.tilman_neumann.jml.factor.siqs.tdiv.TDivReport;
import de.tilman_neumann.jml.factor.squfof.SquFoF31;
import de.tilman_neumann.jml.factor.squfof.SquFoF63;
import de.tilman_neumann.jml.primes.probable.BPSWTest;
import de.tilman_neumann.util.SortedMultiset;
import de.tilman_neumann.util.Timer;

/**
 * A trial division engine where partials can have several large factors.
 * Division is carried out using UnsignedBigInt; this way less intermediate objects are created.
 * 
 * @author Tilman Neumann
 */
public class TDiv_LQS_smallq implements TDiv_LQS {
	private static final Logger LOG = LogManager.getLogger(TDiv_LQS_smallq.class);
	private static final boolean DEBUG = false;
	
	// factor argument and polynomial parameters
	private BigInteger kN;
	private int special_q;
	private BQF_xy bqf;
	
	/** Q is sufficiently smooth if the unfactored Q_rest is smaller than this bound depending on N */
	private double smoothBound;

	// prime base
	private int[] primes;
	private int baseSize;
	private int pMax;
	private BigInteger pMaxSquare;

	/** buffers for trial division engine. */
	private UnsignedBigInt Q_rest_UBI = new UnsignedBigInt(new int[50]);
	private UnsignedBigInt quotient_UBI = new UnsignedBigInt(new int[50]);

	private BPSWTest probablePrimeTest;

	private SquFoF31 squFoF31; // used for Q <= 2^52 that pass trial division
	private SquFoF63 squFoF63; // used for 2^53 <= Q <= 2^59
	private SIQS_Small qsInternal; // Nested SIQS for Q_rest >= 2^60.
	
	// result: two arrays that are reused, their content is _copied_ to AQ-pairs
	private SortedIntegerArray smallFactors = new SortedIntegerArray();
	private SortedLongArray bigFactors = new SortedLongArray();
	private AQPairFactory aqPairFactory = new AQPairFactory();

	// statistics
	private Timer timer = new Timer();
	private long testCount, sufficientSmoothCount;
	private long duration;

	public TDiv_LQS_smallq() {
		this.probablePrimeTest = new BPSWTest();
		this.squFoF31 = new SquFoF31();
		this.squFoF63 = new SquFoF63();
		this.qsInternal = new SIQS_Small(0.305F, 0.37F, null, new SIQSPolyGenerator(), 10, false);
	}

	public String getName() {
		return "TDiv_smallq";
	}

	public void initializeForN(int k, BigInteger N, double N_dbl, BigInteger kN, double smoothBound, int[] primesArray, int baseSize) {
		this.kN = kN;
		// the biggest unfactored rest where some Q is considered smooth enough for a congruence.
		this.smoothBound = smoothBound;
		if (DEBUG) LOG.debug("smoothBound = " + smoothBound + " (" + (64 - Long.numberOfLeadingZeros((long)smoothBound)) + " bits)");
		
		// prime base
		this.primes = primesArray;
		this.baseSize = baseSize;
		pMax = primes[baseSize-1];
		pMaxSquare = BigInteger.valueOf(pMax * (long) pMax);

		// statistics
		if (ANALYZE) testCount = sufficientSmoothCount = 0;
		if (ANALYZE) this.duration = 0;
	}

	public void initializeForSpecialQ(int special_q, BQF_xy bqf) {
		this.special_q = special_q;
		this.bqf = bqf;
	}

	public List<AQPair> testList(ArrayList<IntPair> xyList) {
		if (ANALYZE) timer.capture();

		// do trial division with sieve result
		ArrayList<AQPair> aqPairs = new ArrayList<AQPair>();
		for (IntPair xy : xyList) {
			smallFactors.reset();
			bigFactors.reset();
			if (ANALYZE) testCount++;
			
			int x = xy.x;
			int y = xy.y;
			BigInteger Q = bqf.evaluate(x, y); // smooth candidate
			AQPair aqPair = test(Q, x, y);
			if (aqPair != null) {
				// Q(x) was found sufficiently smooth to be considered a (partial) congruence
				aqPairs.add(aqPair);
				if (ANALYZE) sufficientSmoothCount++;
				if (DEBUG) {
					BigInteger A = aqPair.getA();
					LOG.debug("kN = " + kN + ": Found congruence " + aqPair);
					assertEquals(A.multiply(A).mod(kN), Q.mod(kN));
					// make sure that the product of factors gives Q
					SortedMultiset<Long> allQFactors = aqPair.getAllQFactors();
					BigInteger testProduct = I_1;
					for (Map.Entry<Long, Integer> entry : allQFactors.entrySet()) {
						BigInteger prime = BigInteger.valueOf(entry.getKey());
						int exponent = entry.getValue();
						testProduct = testProduct.multiply(prime.pow(exponent));
					}
					assertEquals(Q, testProduct);
				}
			}
		}
		if (ANALYZE) duration += timer.capture();
		return aqPairs;
	}
	
	private AQPair test(BigInteger Q, int x, int y) {
		// sign
		BigInteger Q_rest = Q;
		if (Q.signum() < 0) {
			smallFactors.add(-1);
			Q_rest = Q.negate();
		}
		
		// Remove multiples of 2
		int lsb = Q_rest.getLowestSetBit();
		if (lsb > 0) {
			smallFactors.add(2, (short)lsb);
			Q_rest = Q_rest.shiftRight(lsb);
		}

		Q_rest_UBI.set(Q_rest);
		// divide by special q first
		while (true) {
			int rem = Q_rest_UBI.divideAndRemainder(special_q, quotient_UBI);
			if (rem > 0) break;
			// Division was exact. "special" q is a big factor, because it can not be found in the prime base.
			// So between different q, full congruences can be found only if two partials have the same q.
			smallFactors.add(special_q);
			UnsignedBigInt tmp = Q_rest_UBI;
			Q_rest_UBI = quotient_UBI;
			quotient_UBI = tmp;
		}
		
		// TODO: is it possible in the lattice sieve to compare x,y with polynomial roots (mod p) to speed up tdiv?

		for (int pIndex = baseSize-1; pIndex > 0; pIndex--) { // p[0]=2 was already tested
			int p = primes[pIndex];
			while (true) {
				int rem = Q_rest_UBI.divideAndRemainder(p, quotient_UBI);
				if (rem>0) break;
				// remainder == 0 -> the division was exact. assign quotient to Q_rest and add p to factors
				UnsignedBigInt tmp = Q_rest_UBI;
				Q_rest_UBI = quotient_UBI;
				quotient_UBI = tmp;
				smallFactors.add(p);
				if (DEBUG) {
					BigInteger pBig = BigInteger.valueOf(p);
					BigInteger[] div = Q_rest.divideAndRemainder(pBig);
					assertEquals(div[1].intValue(), rem);
					Q_rest = div[0];
				}
			}
		}
		if (Q_rest_UBI.isOne()) {
			// We need the AQ-pair factory, because special_q may be contained more than once.
			if (DEBUG) LOG.debug("Found almost smooth, the only large factor being special q! bigFactors = " + bigFactors);
			BigInteger A = computeCongruentQuadratic(Q, x, y);
			return aqPairFactory.create(A, smallFactors, bigFactors);
		}
		Q_rest = Q_rest_UBI.toBigInteger();
		
		// Division by all p<=pMax was not sufficient to factor Q completely.
		// The remaining Q_rest is either a prime > pMax, or a composite > pMax^2.
		if (Q_rest.doubleValue() >= smoothBound) return null; // Q is not sufficiently smooth
		
		// now we consider Q as sufficiently smooth. then we want to know all prime factors, as long as we do not find one that is too big to be useful.
		if (DEBUG) LOG.debug("test(): pMax=" + pMax + " < Q_rest=" + Q_rest + " < smoothBound=" + smoothBound + " -> resolve all factors");
		boolean isSmooth = factor_recurrent(Q_rest);
		if (DEBUG) if (bigFactors.size()>2) LOG.debug("Found " + bigFactors.size() + " distinct big factors!");
		if (!isSmooth) return null;
		BigInteger A = computeCongruentQuadratic(Q, x, y);
		return aqPairFactory.create(A, smallFactors, bigFactors);
	}

	/**
	 * Compute X such that X^2 == Q(x,y) (mod kN).
	 * In practice this is solved by X = (Ax + By)/g (mod kN) = [(Ax + By) * (g^-1 mod kN)] mod kN.
	 * 
	 * @param Q
	 * @param x
	 * @param y
	 * @return X such that X^2 == Q(x,y) (mod kN)
	 */
	BigInteger computeCongruentQuadratic(BigInteger Q, int x, int y) {
		BigInteger Ax = bqf.A.multiply(BigInteger.valueOf(x));
		BigInteger By = bqf.B.multiply(BigInteger.valueOf(y));
		BigInteger X = Ax.add(By).multiply(bqf.gInv).mod(kN); // precomputed g^-1 (mod kN) makes this quite fast
		if (DEBUG) LOG.info("Q=" + Q + ", x=" + x + ", y = " + y + ", bqf = " + bqf + ", g = " + bqf.g + ", gInv = " + bqf.gInv + " -> X = " + X);
		return X;
	}
	
	private boolean factor_recurrent(BigInteger Q_rest) {
		if (Q_rest.compareTo(pMaxSquare)<0) {
			// we divided Q_rest by all primes <= pMax and we have Q_rest < pMax^2 -> it must be prime
			if (DEBUG) {
				LOG.debug("factor_recurrent(): Q_rest = " + Q_rest + " < pMax^2 -> Q_rest is prime");
				assertTrue(probablePrimeTest.isProbablePrime(Q_rest));
			}
			if (Q_rest.bitLength() > 31) return false;
			bigFactors.add(Q_rest.intValue());
			return true;
		}
		// now we can not do without isProbablePrime(), because calling findSingleFactor() may not return when called with a prime argument
		if (probablePrimeTest.isProbablePrime(Q_rest)) {
			// Q_rest is a (probable) prime >= pMax^2. Such big factors do not help to find smooth congruences, so we ignore the partial.
			if (DEBUG) LOG.debug("factor_recurrent(): Q_rest = " + Q_rest + " is probable prime > pMax^2 -> ignore");
			return false;
		} // else: Q_rest is surely not prime
		
		// Find a factor of Q_rest, where Q_rest is odd and has two+ factors, each greater than pMax.
		// At N with 200 bit we have pMax ~ 17 bit, thus Q_rest >= 34 bit -> trial division is no help here.
		BigInteger factor1;
		int Q_rest_bits = Q_rest.bitLength();
		if (Q_rest_bits < 53) {
			if (DEBUG) LOG.debug("factor_recurrent(): pMax^2 = " + pMaxSquare + ", Q_rest = " + Q_rest + " (" + Q_rest_bits + " bits) not prime -> use squFoF31");
			factor1 = squFoF31.findSingleFactor(Q_rest);
		} else if (Q_rest_bits < 60) {
			if (DEBUG) LOG.debug("factor_recurrent(): pMax^2 = " + pMaxSquare + ", Q_rest = " + Q_rest + " (" + Q_rest_bits + " bits) not prime -> use squFoF63");
			factor1 = squFoF63.findSingleFactor(Q_rest);
		} else {
			if (DEBUG) LOG.debug("factor_recurrent(): pMax^2 = " + pMaxSquare + ", Q_rest = " + Q_rest + " (" + Q_rest_bits + " bits) not prime -> use qsInternal");
			factor1 = qsInternal.findSingleFactor(Q_rest);
		}
		// Here we can not exclude factors > 31 bit because they may have 2 prime factors themselves.
		BigInteger factor2 = Q_rest.divide(factor1);
		if (DEBUG) LOG.debug("factor_recurrent(): Q_rest = " + Q_rest + " (" + Q_rest_bits + " bits) = " + factor1 + " * " + factor2);
		return factor_recurrent(factor1) && factor_recurrent(factor2);
	}

	public TDivReport getReport() {
		return new TDivReport(testCount, sufficientSmoothCount, duration, 0, 0, 0, 0);
	}
	
	public void cleanUp() {
		primes = null;
		qsInternal.cleanUp();
	}
}
