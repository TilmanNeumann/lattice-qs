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
import de.tilman_neumann.jml.factor.siqs.tdiv.TDivReport;
import de.tilman_neumann.util.SortedMultiset;
import de.tilman_neumann.util.Timer;

/**
 * A trial division engine where partials can have several large factors.
 * Division is carried out using UnsignedBigInt; this way less intermediate objects are created.
 * 
 * @author Tilman Neumann
 */
public class TDiv_LQS_1L_bigq implements TDiv_LQS {
	private static final Logger LOG = LogManager.getLogger(TDiv_LQS_1L_bigq.class);
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

	/** buffers for trial division engine. */
	private UnsignedBigInt Q_rest_UBI = new UnsignedBigInt(new int[50]);
	private UnsignedBigInt quotient_UBI = new UnsignedBigInt(new int[50]);

	// result: two arrays that are reused, their content is _copied_ to AQ-pairs
	private SortedIntegerArray smallFactors = new SortedIntegerArray();
	private SortedLongArray bigFactors = new SortedLongArray();
	private AQPairFactory aqPairFactory = new AQPairFactory();

	// statistics
	private Timer timer = new Timer();
	private long testCount, sufficientSmoothCount;
	private long duration;

	public TDiv_LQS_1L_bigq() {
	}

	/* (non-Javadoc)
	 * @see de.tilman_neumann.lqs.TDiv_LQS#getName()
	 */
	@Override
	public String getName() {
		return "TDiv_1L_bigq";
	}

	/* (non-Javadoc)
	 * @see de.tilman_neumann.lqs.TDiv_LQS#initializeForN(int, java.math.BigInteger, double, java.math.BigInteger, double, int[], int, boolean)
	 */
	@Override
	public void initializeForN(int k, BigInteger N, double N_dbl, BigInteger kN, double smoothBound, int[] primesArray, int baseSize) {
		this.kN = kN;
		// the biggest unfactored rest where some Q is considered smooth enough for a congruence.
		this.smoothBound = smoothBound;
		if (DEBUG) LOG.debug("smoothBound = " + smoothBound + " (" + (64 - Long.numberOfLeadingZeros((long)smoothBound)) + " bits)");
		
		// prime base
		this.primes = primesArray;
		this.baseSize = baseSize;

		// statistics
		if (ANALYZE) testCount = sufficientSmoothCount = 0;
		if (ANALYZE) this.duration = 0;
	}

	/* (non-Javadoc)
	 * @see de.tilman_neumann.lqs.TDiv_LQS#initializeForSpecialQ(int, de.tilman_neumann.lqs.BQF_xy)
	 */
	@Override
	public void initializeForSpecialQ(int special_q, BQF_xy bqf) {
		this.special_q = special_q;
		this.bqf = bqf;
	}

	/* (non-Javadoc)
	 * @see de.tilman_neumann.lqs.TDiv_LQS#testList(java.util.ArrayList)
	 */
	@Override
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
			bigFactors.add(special_q);
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
		if (Q_rest.bitLength()>31 || Q_rest.doubleValue() >= smoothBound) return null; // Q is not sufficiently smooth
		
		// Q is sufficiently smooth
		if (DEBUG) LOG.debug("Sufficient smooth big factor = " + Q_rest);
		bigFactors.add(Q_rest.intValue());
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

	/* (non-Javadoc)
	 * @see de.tilman_neumann.lqs.TDiv_LQS#getReport()
	 */
	@Override
	public TDivReport getReport() {
		return new TDivReport(testCount, sufficientSmoothCount, duration, 0, 0, 0, 0);
	}
	
	/* (non-Javadoc)
	 * @see de.tilman_neumann.lqs.TDiv_LQS#cleanUp()
	 */
	@Override
	public void cleanUp() {
		primes = null;
	}
}
