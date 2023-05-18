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

import org.apache.log4j.Logger;

import de.tilman_neumann.jml.base.UnsignedBigInt;
import de.tilman_neumann.jml.gcd.EEA31;

import static org.junit.Assert.*;

/**
 * A binary quadratic form Q(x, y) = A*x^2 + 2Bxy + C*y^2 having discriminant D = 4*(B^2 - AC).
 * 
 * @author Tilman Neumann
 */
public class BQF {
	private static final Logger LOG = Logger.getLogger(BQF.class);
	private static final boolean DEBUG = false;
	
	private EEA31 eea31 = new EEA31();
	
	// factor arguments
	int k;
	BigInteger N, kN, kN4;
	
	// the coefficients and variables derived from them
	public BigInteger A, B, twoB, C;
	private UnsignedBigInt A_UBI, B_UBI;
	
	/**
	 * Full constructor.
	 * @param k Knuth-Schroeppel multiplier
	 * @param N the original factoring target
	 * @param kN k*N
	 * @param g the sqrt of A (Pomerance 85 name convention)
	 * @param A
	 * @param B
	 * @param C
	 */
	public BQF(int k, BigInteger N, BigInteger kN, BigInteger A, BigInteger B, BigInteger C) {
		this.k = k;
		this.N = N;
		this.kN = kN;
		this.kN4 = kN.shiftLeft(2);
		this.A = A;
		this.B = B;
		this.twoB = B.shiftLeft(1);
		this.C = C;
		// In the secondary BQF Q'(e, f) = a*e^2 + 2bef + c*f^2, all coefficients a, b, c can be negative.
		// Compute UnsignedBigInts only for positive coefficients:
		this.A_UBI = A.signum()>0 ? new UnsignedBigInt(A) : null;
		this.B_UBI = B.signum()>0 ? new UnsignedBigInt(B) : null;
		
		if (DEBUG) {
			// test discriminant/4
			assertEquals(kN, B.multiply(B).subtract(A.multiply(C)));
		}
	}
	
	/**
	 * @return discriminant D = 4 (B^2-AC) = 4kN
	 */
	public BigInteger discriminant() {
		return kN4;
	}
	
	/**
	 * Evaluate this BQF at (x, y)
	 * @param x
	 * @param y
	 * @return Q(x, y)
	 */
	public BigInteger evaluate(int x, int y) {
		BigInteger xBigSquare = BigInteger.valueOf(x *(long)x);
		BigInteger xyBig = BigInteger.valueOf(x *(long)y);
		BigInteger yBigSquare = BigInteger.valueOf(y *(long)y);
		BigInteger ATerm = A.multiply(xBigSquare);
		BigInteger BTerm = twoB.multiply(xyBig);
		BigInteger CTerm = C.multiply(yBigSquare);
		return ATerm.add(BTerm).add(CTerm);
	}
	
	/**
	 * Solve Q(x,y) == 0 (mod p).
	 * @param p
	 * @param t modular sqrt of kN (mod p)
	 * @return two roots, or null if the modular inverse inverse (1/A) (mod p) does not exist
	 */
	public int[] solvePolynomialRootsWithYEquals1(int p, int t) {
		// We want to solve Ax^2 + 2Bxy + Cy^2 == 0 (mod p) with A, B, C given and y=1.
		// <=> Ax^2 + 2Bx == -C (mod p) // * A // + B^2
		// <=> (Ax)^2 + 2(Ax)B + B^2 == B^2 - AC (mod p)
		// <=> (Ax + B)^2 == B^2 - AC (mod p)
		// <=> (Ax + B)^2 == kN (mod p)
		// The t with t^2 == kN (mod p) is already given...
		
		// Resolve +-t = Ax + B (mod p) towards x:
		// In the primary BQF we have A a large prime -> the modular inverse inverse (1/A) (mod p) should exist always;
		// in the secondary BQF we need to check it.
		BigInteger pBig = BigInteger.valueOf(p);
		int Amodp = A_UBI!=null ? A_UBI.mod(p) : A.mod(pBig).intValue();
		// TODO: In the secondary BQF, A and B can be negative. We can not use UnsignedBigInt then. Can that be cured?
		if (Amodp == 0) {
			// the modular inverse inverse (1/A) (mod p) does not exist
			//LOG.debug("the modular inverse (1/" + A + " (mod " + p + ") does not exist");
			//assertEquals(A.mod(BigInteger.valueOf(p)).intValue(), Amodp);
			return null;
		}
		
		// First solution is x1 = (t-B)/A (mod p)
		long Ainvp = eea31.modularInverse(Amodp, p);
		int Bmodp = B_UBI!=null ? B_UBI.mod(p) : B.mod(pBig).intValue();
		int t_minus_Bmodp = t-Bmodp; // TODO: Here we could save some operations in case t=0
		if (t_minus_Bmodp < 0) t_minus_Bmodp += p;
		int x1 = (int) ((t_minus_Bmodp * Ainvp) % p);
		if (DEBUG) {
			LOG.debug("p = " + p + ", t = " + t + ", A = " + A + ", B = " + B + ", Amodp = " + Amodp + ", Bmodp = " + Bmodp);
			LOG.debug("Ainvp = " + Ainvp + ", t_minus_Bmodp=" + t_minus_Bmodp + ": x1 = " + x1);
//			BigInteger pBig = BigInteger.valueOf(p);
			assertEquals(BigInteger.valueOf(Amodp).modInverse(pBig).intValue(), Ainvp); // wrong if the modular inverse (1/A) mod p does not exist
			assertTrue(0 <= t_minus_Bmodp);
			assertTrue(t_minus_Bmodp < p);
			int kNmodp = kN.mod(pBig).intValue();
			assertEquals(kNmodp, A.multiply(BigInteger.valueOf(x1)).add(B).pow(2).mod(pBig).intValue());
			assertEquals(kNmodp, (t*(long)t) % p);
		}
		if (t==0) return new int[] {x1, x1}; // only one solution
		
		// Compute second solution x2 = (p-t-B)/A (mod p)
		int minus_t_minus_Bmodp = p-t-Bmodp;
		if (minus_t_minus_Bmodp < 0) minus_t_minus_Bmodp += p;
		int x2 = (int) ((minus_t_minus_Bmodp * Ainvp) % p);
		if (DEBUG) {
			LOG.debug("minus_t_minus_Bmodp=" + minus_t_minus_Bmodp + ": x2 = " + x2);
//			BigInteger pBig = BigInteger.valueOf(p);
			assertTrue(0 <= minus_t_minus_Bmodp);
			assertTrue(minus_t_minus_Bmodp < p);
			int kNmodp = kN.mod(pBig).intValue();
			assertEquals(kNmodp, A.multiply(BigInteger.valueOf(x2)).add(B).pow(2).mod(pBig).intValue());
		}
		return new int[] {x1, x2};
	}

	@Override
	public String toString() {
		return "BQF for discriminant " + discriminant() + " = (" + A + ", " + B + ", " + C + ")";
	}
}
