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

import java.math.BigInteger;
import java.util.HashSet;

import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

import de.tilman_neumann.jml.factor.FactorException;
import de.tilman_neumann.jml.modular.JacobiSymbol;
import de.tilman_neumann.jml.modular.ModularSqrt_BB;
import de.tilman_neumann.jml.primes.probable.BPSWTest;
import de.tilman_neumann.jml.roots.Roots;

import static de.tilman_neumann.jml.base.BigIntConstants.*;
import static org.junit.Assert.*;

/**
 * Factory to construct a bivariate polynomial Q(x,y) = Ax^2 + 2Bxy + Cy^2 such that we get quadratic residues from (Ax+By)/g (mod kN).
 * This is very similar to the approach in [Pomerance 1985 "The quadratic sieve factoring algorithm"]
 * 
 * @author Tilman Neumann
 */
public class BQFFactory {
	private static final Logger LOG = LogManager.getLogger(BQFFactory.class);
	private static final boolean DEBUG = false;

	/**
	 * How to compute the b-parameter:
	 * 0 = simple formula doing arithmetics (mod p^2)
	 * 1 = more elaborated and faster version that needs arithmetics (mod p) only
	 */
	private static final int B_COMPUTATION_STYLE = 1;

	private BPSWTest bpsw = new BPSWTest();
	private JacobiSymbol jacobiEngine = new JacobiSymbol();
	private ModularSqrt_BB modEngine = new ModularSqrt_BB();

	/**
	 * Construct a binary quadratic form Q(x, y) = A*x^2 + 2Bxy + C*y2 with discriminant D = 4(B^2-AC) = 4kN.
	 * 
	 * @param k
	 * @param N
	 * @param kN
	 * @return BQF (A,B,C)
	 * @throws FactorException in the rare case where we found some A with Jacobi(4kN, A) == 0
	 */
	// [Crandall/Pomerance p.239 ff.]:
	// BQF with negative discriminant
	// * are easier to handle
	// * do not represent both positive and negative numbers
	// * have A>0 which forces C>0
	// So far I have D>0, A>0, C<0.
	// TODO: Try D<0, i.e. negative discriminants ?
	public BQF_xy createBQFWithDiscriminant4kN(int k, BigInteger N, BigInteger kN, int[] primesArray, int primeBaseSize, long sieveArrayArea) throws FactorException {

		// required to effectively check that g is not in the prime base
		BigInteger maxPrime_big = BigInteger.valueOf(primesArray[primeBaseSize-1]);
		HashSet<Integer> hashedPrimeBase = new HashSet<Integer>(); // TODO is computed again in LatticeVectorSieve
		for (int i=0; i<primeBaseSize; i++) {
			hashedPrimeBase.add(primesArray[i]);
		}

		// find some initial g ~ sqrt(sqrt(2kN)/sieveArrayArea); it does not need to be prime yet.
//		BigInteger sieveArraySize_big = BigInteger.valueOf(sieveArrayArea);
//		BigInteger sieveArraySizeSquare_big = sieveArraySize_big.multiply(sieveArraySize_big);
//		BigInteger g = Roots.ithRoot(kN.shiftLeft(1).divide(sieveArraySizeSquare_big), 4)[0]; // Pomerance's proposal for MPQS
		BigInteger g = Roots.ithRoot(kN, 4)[0]; // size-balanced
		// TODO optimal size of g? size-balanced choice was faster at first glance
		
		// Find a prime g such that
		// a) g == 3 (mod 4), not required here but good for performance
		// b) kN is a quadratic residue (mod g)
		// c) g is not in the prime base
		if (DEBUG) LOG.debug("pMax = " + maxPrime_big + ", start searching for g at g=" + g);
		while (true) {
			// [Silverman 1987: "The multiple polynomial quadratic sieve", p.332]:
			// "It is sufficient for practical purposes that D be only a probable prime" (Silvermans D = our g)
			g = bpsw.nextProbablePrime(g);
			
			// it is not required to choose g == 3 (mod 4) but very good for performance
			if (g.and(I_3).intValue()!=3) {
				if (DEBUG) LOG.debug("g = " + g + " != 3 (mod 4) -> skip");
				continue;
			}
			
			if (g.compareTo(maxPrime_big)<=0 && hashedPrimeBase.contains(g.intValue())) {
				if (DEBUG) LOG.debug("g = " + g + " is contained in prime base -> skip");
				continue;
			}
			// otherwise g looks ok so far -> check if kN is a quadratic residue (mod g)
			int jacobi = jacobiEngine.jacobiSymbol(kN, g);
			if (jacobi>0) {
				// kN is quadratic residue (mod g) -> that's what we were looking for.
				if (DEBUG) LOG.debug("found g = " + g);
				break; // use this g
			}
			if (jacobi==0) {
				// g divides kN -> Since g is typically much bigger than k, this means we found a factor of N! -> test gcd
				if (DEBUG) LOG.debug("g = " + g + " divides kN -> a factor of N?");
				BigInteger gcd = N.gcd(g);
				if (gcd.compareTo(I_1)>0 && gcd.compareTo(N)<0) throw new FactorException(gcd);
			} else {
				if (DEBUG) LOG.debug("kN = " + k + " * " + N + " is not a quadratic residue (mod g), g = " + g + " -> skip");
			}
		}
		// a = g^2
		BigInteger a = g.multiply(g);
		
		// Compute b such that b^2 == kN (mod a). this requires g == 3 (mod 4) as chosen above.
		BigInteger b;
		switch (B_COMPUTATION_STYLE) {
		case 0:
			b = modEngine.modularSqrtModSquare_v01(kN, g, a);
			break;
		case 1:
			b = modEngine.modularSqrtModSquare_v02(kN, g, a);
			break;
		default: throw new IllegalStateException("Unknonwn b-computation style: " + B_COMPUTATION_STYLE);
		}
		if (DEBUG) assertEquals(I_0, b.multiply(b).subtract(kN).mod(a));
		
		BigInteger c = b.multiply(b).subtract(kN).divide(a);
		BQF_xy bqf = new BQF_xy(k, N, kN, g, a, b, c);
		if (DEBUG) assertEquals(kN.shiftLeft(2), bqf.discriminant());
		return bqf;
	}
}
