/*
 * lattics-qs (LQS) is a study of a lattice quadratic sieve proposed by R.D.Silverman in https://www.mersenneforum.org/showthread.php?t=14080.
 * Copyright (C) 2018 Tilman Neumann (www.tilman-neumann.de)
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

import static de.tilman_neumann.jml.base.BigIntConstants.ZERO;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashSet;

import org.apache.log4j.Logger;

import de.tilman_neumann.jml.base.UnsignedBigInt;
import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.factor.FactorException;
import de.tilman_neumann.jml.factor.siqs.sieve.SieveReport;
import de.tilman_neumann.jml.modular.ModularSqrt31;
import de.tilman_neumann.util.SortedMultiset;
import de.tilman_neumann.util.Timer;

/**
 * A lattice vector sieve (with special-q), sieving in the (e, f) plane.
 * 
 * This implementation may be very bad or not even what it is supposed to be, because I realized it without reading
 * [Pollard 1993: "The lattice sieve"]. I could not find it for free in the internet and was too lazy yet to go to a library.
 * My main source of information was [Franke/Kleinjung 2005: "Continued fractions and lattice sieving"]
 * 
 * @author Tilman Neumann
 */
// TODO implement segmented sieve
// TODO implement line sieve
// TODO For large numbers, we might need to convert (int*int) multiplications to (int*long) multiplications. At 210 bit, this had no impact yet, though...
public class LatticeVectorSieve implements LatticeSieve {
	private static final Logger LOG = Logger.getLogger(LatticeVectorSieve.class);
	private static final boolean DEBUG = false;
	private static final boolean TEST_FACTORS = false;
	private static final boolean TEST_SIEVE_ARRAY = false;

	// factor arguments
	int k;
	BigInteger N, kN;
	UnsignedBigInt kN_UBI;
	
	// prime base
	int[] primesArray;
	HashSet<Integer> hashedPrimeBase;
	int primeBaseSize;
	int pMax;
	int pMinIndex;
	
	// polynomial
	private ModularSqrt31 modularSqrtEngine = new ModularSqrt31();
	private int[] tArray;
	private BQF Qxy;

	// profiling
	private boolean profile;
	private Timer timer = new Timer();
	private long initDuration, sieveDuration, collectDuration;

	@Override
	public String getName() {
		return "vectorSieve";
	}
	
	@Override
	public void initializeForN(int k, BigInteger N, BigInteger kN, int[] primesArray, int[] tArray, int primeBaseSize, LQSSieveParams sieveParams, BQF_xy Qxy, boolean profile) {
		this.k = k;
		this.N = N;
		this.kN = kN;
		this.kN_UBI = new UnsignedBigInt(kN);
		this.primesArray = primesArray;
		this.tArray = tArray;
		this.primeBaseSize = primeBaseSize;
		this.pMinIndex = sieveParams.pMinIndex;
		
		// just for debugging
		hashedPrimeBase = new HashSet<Integer>();
		for (int p : primesArray) {
			hashedPrimeBase.add(p);
		}
		pMax = primesArray[primeBaseSize-1];
		
		this.Qxy = Qxy;

		// profiling
		this.profile = profile;
		initDuration = sieveDuration = collectDuration = 0;
	}

	@Override
	public ArrayList<IntPair> sieve(int q, byte[] logPArray, byte[][] sieveArray, byte[] initializedSieveLine, byte[][] dontUseArray, int sieveArraySideLength) throws FactorException {
		if (profile) timer.capture();

		// Compute modular sqrt of kN (mod q)
		int kN_mod_q = kN_UBI.mod(q);
		int t = (kN_mod_q > 0) ? modularSqrtEngine.modularSqrt(kN_mod_q, q) : 0;
		
		// Fix y (e.g. to y=1) and solve Q(x, y) == 0 (mod q). This gives one or two roots.
		int[] qRoots = Qxy.solvePolynomialRootsWithYEquals1(q, t);
		if (DEBUG) {
			LOG.debug("special q = " + q + ", qRoot1 = " + qRoots[0] + ", qRoot2 = " + qRoots[1]);
			assertTrue(qRoots!=null);
			assertTrue(qRoots[0] != qRoots[1]);
			// The roots must satisfy Q(root, 1) == 0 (mod q)
			BigInteger qBig = BigInteger.valueOf(q);
			BigInteger Q1 = Qxy.evaluate(qRoots[0], 1);
			assertTrue(Q1.mod(qBig).equals(ZERO));
			BigInteger Q2 = Qxy.evaluate(qRoots[1], 1);
			assertTrue(Q2.mod(qBig).equals(ZERO));
		}
		
		// Silverman: "Find reduced lattices for [q r1][0 1] and [q r2][0 1]. This yields two sets of vectors V1, V2 and W1, W2."
		// An unreduced lattice is given by linear combinations of the first solution vector (r, 1) and vector (q, 0).
		LatticeBase qBase1 = new LatticeBase(qRoots[0], 1, q, 0).reduce();
		if (DEBUG) {
			LOG.debug("q=" + q + ", qRoot1=" + qRoots[0] + ": qBase1 = " + qBase1 + " has det = " + qBase1.determinant());
			// Test reduced base for e=2, f=2
			int x = 2*qBase1.b10 + 2*qBase1.b20;
			int y = 2*qBase1.b11 + 2*qBase1.b21;
			BigInteger Q = Qxy.evaluate(x, y);
			assertTrue(Q.mod(BigInteger.valueOf(q)).intValue() == 0);
		}

		// sieve loop for q-reduced base qBase1
		ArrayList<IntPair> latticeSieveSmoothXYList = sieveOneQRoot(q, qBase1, logPArray, sieveArray, initializedSieveLine, dontUseArray, sieveArraySideLength);

		if (qRoots[1] != qRoots[0]) {
			LatticeBase qBase2 = new LatticeBase(qRoots[1], 1, q, 0).reduce();
			if (DEBUG) {
				LOG.debug("q=" + q + ", qRoot2=" + qRoots[1] + ": qBase2 = " + qBase2 + " has det = " + qBase2.determinant());
				// Test reduced base for e=2, f=2
				int x = 2*qBase2.b10 + 2*qBase2.b20;
				int y = 2*qBase2.b11 + 2*qBase2.b21;
				BigInteger Q = Qxy.evaluate(x, y);
				assertTrue(Q.mod(BigInteger.valueOf(q)).intValue() == 0);
			}
	
			// sieve loop for q-reduced base qBase2
			ArrayList<IntPair> latticeSieveSmoothXYList2 = sieveOneQRoot(q, qBase2, logPArray, sieveArray, initializedSieveLine, dontUseArray, sieveArraySideLength);
			latticeSieveSmoothXYList.addAll(latticeSieveSmoothXYList2);
		}
		
		return latticeSieveSmoothXYList;
	}
	
	/**
	 * Do lattice sieving with one of the two q-reduced lattice bases.
	 * @param q the "special q"
	 * @param qBase reduced base of the lattice generated by q, or more exactly by Q(x, y) == 0 (mod q)
	 * @param logPArray
	 * @param sieveArray a quadratic 2D sieve array having side length sieveArraySideLength
	 * @param initializedSieveLine
	 * @param dontUseArray
	 * @param sieveArraySideLength
	 * @return list of (x,y)-pairs representing smooth candidates
	 */
	private ArrayList<IntPair> sieveOneQRoot(
			int q, LatticeBase qBase, byte[] logPArray, byte[][] sieveArray, byte[] initializedSieveLine, byte[][] dontUseArray, int sieveArraySideLength) throws FactorException {
		
		// compute second BQF in (e, f); here we may get coefficients a, b < 0
		BQF Qef = computeSecondBQF(q, qBase);
		
		// Use a quadratic sieve array as Silverman did.
		// Allowing both x and y to attain negative values would be no help, because then the polynomial would produce the same values twice.
		// Thus we restrict y-values to be positive.
		int halfSieveArraySideLength = sieveArraySideLength>>1;
		int eMin = -halfSieveArraySideLength;
		int eMax = halfSieveArraySideLength-1;
		int fMin = 0;
		int fMax = sieveArraySideLength-1;
		
		// (re-)initialize sieveArray
		for (int i=sieveArraySideLength-1; i>=0; i--) {
			System.arraycopy(initializedSieveLine, 0, sieveArray[i], 0, sieveArraySideLength);
		}
		
		// For each p, construct a reduced basis of 2-dim. vectors for the sublattice of L generated by (pq).
		// There may be one or two distinct lattices for each p.
		// Not sieving with small primes speeds up sieving by about 20% for large N.
		for (int i=primeBaseSize-1; i>=pMinIndex; i--) {
			// Find polynomial roots of Q'(e,f) == 0 (mod p).
			// If the modular inverse (1/A) mod p does not exist, then usually there is no solution to Q'(e,f) == 0 (mod p).
			// In that case here we would erroneously get roots [0, 0], followed by Q1, Q2 (mod p) != 0 below...
			int p = primesArray[i];
			int[] qpRoots = Qef.solvePolynomialRootsWithYEquals1(p, tArray[i]);
			if (qpRoots == null) continue;
			if (DEBUG) {
				LOG.debug("qpRoots for p = " + p + ", t = " + tArray[i] + ": " + qpRoots[0] + ", " + qpRoots[1]);
				// test qpRoots
				BigInteger pBig = BigInteger.valueOf(p);
				BigInteger Q1 = Qef.evaluate(qpRoots[0], 1);
				int Q1modp = Q1.mod(pBig).intValue();
				LOG.debug("q=" + q + ", p=" + p + ", Q1=" + Q1 + ", Q1modp=" + Q1modp);
				assertEquals(0, Q1modp);
				BigInteger Q2 = Qef.evaluate(qpRoots[1], 1);
				int Q2modp = Q2.mod(pBig).intValue();
				LOG.debug("q=" + q + ", p=" + p + ", Q2=" + Q2 + ", Q2modp=" + Q2modp);
				assertEquals(0, Q2modp);
			}
			
			LatticeBase qpBase1 = new LatticeBase(qpRoots[0], 1, p, 0).reduce();
			if (DEBUG) {
				LOG.debug("qpBase1 = " + qpBase1);
				// Test reduced base for c=2, d=2
				int e = 2*qpBase1.b10 + 2*qpBase1.b20;
				int f = 2*qpBase1.b11 + 2*qpBase1.b21;
				int x = e*qBase.b10 + f*qBase.b20;
				int y = e*qBase.b11 + f*qBase.b21;
				BigInteger Q = Qxy.evaluate(x, y);
				assertTrue(Q.mod(BigInteger.valueOf(q)).intValue() == 0);
				assertTrue(Q.mod(BigInteger.valueOf(p)).intValue() == 0); // works for exactly orthogonal sublattice case, too
				//assertTrue(qpBase1.b10 != 0); // wrong, we still get exactly orthogonal sublattices
			}
			if (profile) initDuration += timer.capture();
			
			byte logP = logPArray[i];
			sieveOneP(qBase, p, qpBase1, eMin, eMax, fMin, fMax, sieveArray, sieveArraySideLength, halfSieveArraySideLength, logP);
			if (profile) sieveDuration += timer.capture();

			if (qpRoots[0] != qpRoots[1]) {
				// There is a second solution
				LatticeBase qpBase2 = new LatticeBase(qpRoots[1], 1, p, 0).reduce();
				if (DEBUG) {
					LOG.debug("qpBase2 = " + qpBase2);
					// Test reduced base for c=2, d=2
					int e = 2*qpBase2.b10 + 2*qpBase2.b20;
					int f = 2*qpBase2.b11 + 2*qpBase2.b21;
					int x = e*qBase.b10 + f*qBase.b20;
					int y = e*qBase.b11 + f*qBase.b21;
					BigInteger Q = Qxy.evaluate(x, y);
					assertTrue(Q.mod(BigInteger.valueOf(q)).intValue() == 0);
					assertTrue(Q.mod(BigInteger.valueOf(p)).intValue() == 0); // works for exactly orthogonal sublattice case, too
					//assertTrue(qpBase2.b10 != 0); // wrong, we still get exactly orthogonal sublattices
				}
				if (profile) initDuration += timer.capture();
				
				sieveOneP(qBase, p, qpBase2, eMin, eMax, fMin, fMax, sieveArray, sieveArraySideLength, halfSieveArraySideLength, logP);
				if (profile) sieveDuration += timer.capture();
			} else {
				if (DEBUG) LOG.info("p = " + p + " has only one sublattice!");
			}
		}
		
		// Collect sieve results
		ArrayList<IntPair> smoothCandidates = collectOneRoot(q, qBase, eMin, eMax, fMin, fMax, sieveArray, dontUseArray, halfSieveArraySideLength);
		if (profile) collectDuration += timer.capture();
		return smoothCandidates;
	}

	/**
	 * Compute BQF (a, b, c) = Q'(e, f) = Q(e*v1+f*w1, e*v1+f*w2)/q
	 *
	 * @param q special q
	 * @param qBase q-reduced lattice base
	 * @return
	 */
	private BQF computeSecondBQF(int q, LatticeBase qBase) {
		// compute a, b, c params
		BigInteger a = Qxy.evaluate(qBase.b10, qBase.b11);
		BigInteger b1 = Qxy.A.multiply(BigInteger.valueOf(qBase.b10 * (long) qBase.b20));
		BigInteger b2 = Qxy.B.multiply(BigInteger.valueOf(qBase.b10 * (long) qBase.b21));
		BigInteger b3 = Qxy.B.multiply(BigInteger.valueOf(qBase.b11 * (long) qBase.b20));
		BigInteger b4 = Qxy.C.multiply(BigInteger.valueOf(qBase.b11 * (long) qBase.b21));
		BigInteger b = b1.add(b2).add(b3).add(b4);
		BigInteger c = Qxy.evaluate(qBase.b20, qBase.b21);
		// a, b, c are divisible by q
		BigInteger qBig = BigInteger.valueOf(q);
		BigInteger aDivq = a.divide(qBig);
		BigInteger bDivq = b.divide(qBig);
		BigInteger cDivq = c.divide(qBig);
		
		BQF Qef = new BQF(k, N, kN, aDivq, bDivq, cDivq);
		if (DEBUG) {
			LOG.debug("a = " + a + ", b = " + b + ", c = " + c);
			LOG.debug("a/q = " + aDivq + ", b/c = " + bDivq + ", c/q = " + cDivq);
			LOG.debug("A = " + Qxy.A + ", B = " + Qxy.B + ", C = " + Qxy.C);
			assertEquals(a, aDivq.multiply(qBig));
			assertEquals(b, bDivq.multiply(qBig));
			assertEquals(c, cDivq.multiply(qBig));
			assertEquals(Qxy.discriminant(), Qef.discriminant()); // the discriminant is the same as for BQF (A,B,C)
		}
		return Qef;
	}

	/**
	 * Lattice sieve a single prime p.
	 * @param qBase
	 * @param p
	 * @param qpBase
	 * @param eMin
	 * @param eMax
	 * @param fMin
	 * @param fMax
	 * @param sieveArray
	 * @param halfSieveArraySideLength
	 * @param logP
	 */
	private void sieveOneP(LatticeBase qBase, int p, LatticeBase qpBase, int eMin, int eMax, int fMin, int fMax, byte[][] sieveArray, int sieveArraySideLength, int halfSieveArraySideLength, byte logP) {
		int v1 = qpBase.b10;
		int v2 = qpBase.b11;
		int w1 = qpBase.b20;
		int w2 = qpBase.b21;
		
		int cMin, cMax, dMin, dMax;
		boolean isExactlyOrthogonal = false;
		
		// Compute c, d bounds from
		// e = c*v1 + d*w1    =>    eMin <= c*v1 + d*w1 <= eMax
		// f = c*v2 + d*w2    =>    fMin <= c*v2 + d*w2 <= fMax
		// We do not want negative d, because
		// a) if (c, d) gives (e, f) then (-c, -d) gives (-e, -f)
		// b) if (-c, d) gives (e, f) then (c, -d) gives (-e, -f)
		if (v1==0) {
			// if v1 = w2 = 0, then we have a base that is exactly orthogonal.
			// To be a base at all we must have v2, w1 != 0; actually is is always v2, w1 > 0.
			if (DEBUG) {
				LOG.info("qp-reduced lattice base " + qpBase + " for p=" + p + " is exactly orthogonal!");
				assertEquals(0, w2);
				assertTrue(v2 > 0);
				assertTrue(w1 > 0);
			}
			isExactlyOrthogonal = true;
			// We have v1=0, v2>0, w1>0, w2=0.
			// 1. eMin <= d*w1
			// From w1>0: eMin/w1 <= d		<=> dMin = eMin/w1
			// We have dMin<0 always because eMin<0. Since we do not want negative d, we just choose dMin=0.
			dMin = 0;
			// 2. d*w1 <= eMax
			// From w1>0: d <= eMax/w1		<=> dMax = eMax/w1
			dMax = (int) Math.floor(eMax/(float)w1);
			// 3. fMin <= c*v2
			// From v2>0: fMin/v2 <= c 		<=> cMin = fMin/v2
			// Since fMin=0 we have cMin=0.
			cMin = 0;
			// 4. c*v2 <= fMax
			// From v2>0: c <= fMax/v2		<=> cMax = fMax/v2
			cMax = (int) Math.floor(fMax/(float)v2);
			if (DEBUG) {
				LOG.debug("cMin = " + cMin + ", cMax = " + cMax + ", dMin=" + dMin + ", dMax=" + dMax);
				assertTrue(Math.ceil(eMin/(float)w1) <= 0); // original dMin derivation
				assertTrue(dMax >= 0);
			}
		} else {
			cMin = cMax = 0; // dummy, correct values will be computed for each d
			dMin = 0;
			// The distance of "vector lines" is given by the height of the parallelogram spanned by the qp-reduced bases.
			// The number of "vector lines" to sieve should therefore roughly be sieveArraySideLength / distance
			float dist = (float) (Math.abs(v1*w2 - v2*w1)/Math.sqrt(v1*v1 + v2*v2));
			dMax = (int) (sieveArraySideLength / dist);
		}
		
		for (int d=dMin; d<=dMax; d++) {
			if (!isExactlyOrthogonal) {
				// Compute c bounds from given d and
				// e = c*v1 + d*w1    =>    eMin <= c*v1 + d*w1 <= eMax
				// f = c*v2 + d*w2    =>    fMin <= c*v2 + d*w2 <= fMax
				
				double cMin1, cMax1;
				if (v1 < 0) {
					// We have v1<0, v2>=0, w2>=0.
					// 1. eMin <= c*v1 + d*w1
					// <=> eMin - d*w1 <= c*v1
					// From v1<0: (eMin - d*w1)/v1 >= c		<=> cMax = (eMin - d*w1)/v1
					// 2. c*v1 + d*w1 <= eMax
					// <=> c*v1 <= eMax - d*w1
					// From v1<0: c >= (eMax - d*w1)/v1		<=> cMin = (eMax - d*w1)/v1
					// 3. fMin <= c*v2 + d*w2
					// <=> fMin - d*w2 <= c*v2
					// From v2>=0: (fMin - d*w2)/v2 <= e	<=> cMin = (fMin - d*w2)/v2
					// 4. c*v2 + d*w2 <= fMax
					// <=> c*v2 <= fMax - d*w2
					// From v2>=0: c <= (fMax - d*w2)/v2	<=> cMax = (fMax - d*w2)/v2
					cMin1 = (eMax - d*w1)/(double)v1;
					cMax1 = (eMin - d*w1)/(double)v1;
				} else {
					// We have v1>0, v2>=0, w2>=0.
					// 1. eMin <= c*v1 + d*w1
					// <=> eMin - d*w1 <= c*v1
					// From v1>0: (eMin - d*w1)/v1 <= c		<=> cMin = (eMin - d*w1)/v1
					// 2. c*v1 + d*w1 <= eMax
					// <=> c*v1 <= eMax - d*w1
					// From v1>0: c <= (eMax - d*w1)/v1		<=> cMax = (eMax - d*w1)/v1
					// 3. fMin <= c*v2 + d*w2
					// <=> fMin - d*w2 <= c*v2
					// From v2>=0: (fMin - d*w2)/v2 <= c	<=> cMin = (fMin - d*w2)/v2
					// 4. c*v2 + d*w2 <= fMax
					// <=> c*v2 <= fMax - d*w2
					// From v2>=0: c <= (fMax - d*w2)/v2	<=> cMax = (fMax - d*w2)/v2
					cMin1 = (eMin - d*w1)/(double)v1;
					cMax1 = (eMax - d*w1)/(double)v1;
				}
				final double cMin2 = (fMin - d*w2)/(double)v2;
				final double cMax2 = (fMax - d*w2)/(double)v2;
				
				// Since we are talking about bounds, we have to take the bigger cMin and the smaller cMax.
				// Lower bounds are rounded up, upper bounds are truncated.
				cMin = (int) Math.ceil(Math.max(cMin1, cMin2));
				cMax = (int) Math.floor(Math.min(cMax1, cMax2)); // floor is important for correct bounds!
				if (DEBUG) {
					LOG.debug("cMin1=" + cMin1 + ", cMin2=" + cMin2 + ", cMax1=" + cMax1 + ", cMax2=" + cMax2);
					LOG.debug("v1=" + v1 + ", v2=" + v2 + ", w1=" + w1 + ", w2=" + w2);
					LOG.debug("eMin=" + eMin + ", eMax=" + eMax + ", fMin=" + fMin + ", fMax=" + fMax);
					LOG.debug("d=" + d + ", cMin=" + cMin + ", cMax=" + cMax);
				}
			}
			int cCount = cMax - cMin + 1; // cMin, cMax are inclusive
			
			if (cCount > 0) {
				// Sieve through all vectors on the (e,f)-lattice via the linear combination (e,f) = c*v + d*w, with v=(v1, v2), w=(w1, w2).
				//if (DEBUG) LOG.debug("p=" + p + ", d=" + d + ": cMin = " + cMin + ", cMax = " + cMax + ", cCount = " + cCount);
				int eIndex = cMin*v1 + d*w1 + halfSieveArraySideLength;
				int f = cMin*v2 + d*w2;
				for (int c = cMin; c<=cMax; c++, eIndex += v1, f += v2) { // we want to log the correct c-values
					try {
						// Without double precision in the bounds computations above, here we would get an ArrayIndexOutOfBoundsException for N=210 bit.
						sieveArray[eIndex][f] += logP;
						
						if (DEBUG) {
							int e = eIndex - halfSieveArraySideLength;
							//LOG.debug("c=" + c + ", d=" + d + "; e=" + e + ", f=" + f);
							int x = e*qBase.b10 + f*qBase.b20;
							int y = e*qBase.b11 + f*qBase.b21;
							BigInteger Q = Qxy.evaluate(x, y);
							BigInteger pBig = BigInteger.valueOf(p);
							assertTrue(Q.mod(pBig).equals(ZERO)); // works for exactly orthogonal sublattice case, too

							if (TEST_FACTORS) {
								// verify that all small factors of Q are in the prime base
								SortedMultiset<BigInteger> QFactors = FactorAlgorithm.DEFAULT.factor(Q);
								for (BigInteger factor : QFactors.keySet()) {
									if (factor.bitLength()> 31) break; // too big
									int iFactor = factor.intValue();
									if (1<iFactor && iFactor <= pMax && !hashedPrimeBase.contains(iFactor)) {
										LOG.info("p=" + p + ": Q(x,y) = " + Q + " has factor not in prime base: " + iFactor);
									}
								}
							}
						}
					} catch (ArrayIndexOutOfBoundsException ex) {
						LOG.error("Found exception: " + ex, ex);
						LOG.error("qpReducedLatticeBase = " + qpBase);
						LOG.error("cMin = " + cMin + ", cMax = " + cMax + ", dMin = " + dMin + ", dMax = " + dMax);
						LOG.error("eMin = " + eMin + ", eMax = " + eMax + ", fMin = " + fMin + ", fMax = " + fMax);
						int e = eIndex - halfSieveArraySideLength;
						LOG.error("c = " + c + ", d = " + d + ", e = " + e + ", f = " + f);
						throw ex;
					}
				}
			}
		}
	}
	
	/**
	 * Collect sieve results: All (e,f) of the special_q lattice produce Q(x,y) are divisible by q, so we have to run over all (e,f).
	 *
	 * @param q
	 * @param qBase q-reduced lattice base
	 * @param eMin
	 * @param eMax
	 * @param fMin
	 * @param fMax
	 * @param sieveArray
	 * @param dontUseArray
	 * @param halfSieveArraySideLength
	 * @return smooth (x, y) pair candidates
	 */
	private ArrayList<IntPair> collectOneRoot(int q, LatticeBase qBase, int eMin, int eMax, int fMin, int fMax, byte[][] sieveArray, byte[][] dontUseArray, int halfSieveArraySideLength) {
		// Collect sieve results: Requires to run over the special_q lattice
		int v1 = qBase.b10;
		int v2 = qBase.b11;
		int w1 = qBase.b20;
		int w2 = qBase.b21;
		
		ArrayList<IntPair> allSmooths = new ArrayList<IntPair>();
		ArrayList<IntPair> foundSmooths = new ArrayList<IntPair>();
		for (int e = eMin; e<=eMax; e++) {
			final int eIndex = e+halfSieveArraySideLength;
			byte[] dontUseArrayRow = dontUseArray[eIndex];
			byte[] sieveArrayRow = sieveArray[eIndex];
			for (int f=fMin; f<=fMax; f++) {
				// We only want (e,f) != (0,0) and gcd(e, f) == 1.
				// The gcd-condition reduces the number of required relations.
				// TODO Why, are that duplicates? Probably not exact duplicates, or?
				if (dontUseArrayRow[f]==1) continue;
				if ((sieveArrayRow[f] & 0x80) != 0) {
					// found a smooth candidate!
					int x = e*v1 + f*w1;
					int y = e*v2 + f*w2;
					IntPair xyPair = new IntPair(x, y);
					foundSmooths.add(xyPair);
				}
				if (DEBUG) {
					LOG.debug("e=" + e + ", f=" + f + ": dontUse = " + dontUseArray[e+halfSieveArraySideLength][f]);
					assertTrue(e!=0 || f!=0);
					int x = e*v1 + f*w1;
					int y = e*v2 + f*w2;
					assertTrue(x!=0 || y!=0);
					BigInteger Q = Qxy.evaluate(x, y);
					BigInteger qBig = BigInteger.valueOf(q);
					assertTrue(Q.mod(qBig).equals(ZERO));
					if (TEST_SIEVE_ARRAY) {
						// test if the rest is smooth
						BigInteger QRest = Q.divide(qBig);
						SortedMultiset<BigInteger> factors = FactorAlgorithm.DEFAULT.factor(QRest);
						BigInteger biggestFactor = factors.getBiggestElement();
						if (biggestFactor.bitLength()<32 && biggestFactor.intValue()<=pMax) {
							if (DEBUG) LOG.debug("smooth Q=" + Q + ": QRest=" + QRest + " has factors " + factors + ", biggestFactor=" + biggestFactor + ", pMax=" + pMax);
							allSmooths.add(new IntPair(x, y));
						}
					}
				}
			}
		}
		if (TEST_SIEVE_ARRAY) {
			LOG.info("q=" + q + ": sieve array contains " + allSmooths.size() + " smooth Q: " + allSmooths);
		}
		if (DEBUG) LOG.info("q=" + q + ": vector sieve found " + foundSmooths.size() + " smooth Q: " + foundSmooths);	
		return foundSmooths;
	}
	
	@Override
	public SieveReport getReport() {
		return new SieveReport(initDuration, sieveDuration, collectDuration);
	}
}
