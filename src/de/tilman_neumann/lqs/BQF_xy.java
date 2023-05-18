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

/**
 * The polynomial in x, y: A binary quadratic form Q(x, y) = A*x^2 + 2Bxy + C*y^2 having discriminant D = 4*(B^2 - AC).
 * 
 * @author Tilman Neumann
 */
public class BQF_xy extends BQF {
	public BigInteger g;
	public BigInteger gInv; // (1/g) (mod kN)
	
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
	public BQF_xy(int k, BigInteger N, BigInteger kN, BigInteger g, BigInteger A, BigInteger B, BigInteger C) {
		super(k, N, kN, A, B, C);
		this.g = g;
		gInv = g.modInverse(kN); // for all-BigInteger arguments, I have nothing faster than BigInteger.modInverse()
	}
}
