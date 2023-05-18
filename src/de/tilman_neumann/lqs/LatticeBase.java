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

/**
 * A lattice base of 2-dim. vectors b1, b2.
 * @author Tilman Neumann
 */
// TODO: b10[], b11[], b20[], b21[] arrays would permit better performance than a LatticeBase[] ?
public class LatticeBase {
	public int b10, b11, b20, b21;
	
	/**
	 * Standard constructor.
	 * @param b10
	 * @param b11
	 * @param b20
	 * @param b21
	 */
	public LatticeBase(int b10, int b11, int b20, int b21) {
		this.b10 = b10;
		this.b11 = b11;
		this.b20 = b20;
		this.b21 = b21;
	}
	
	/**
	 * Copy constructor.
	 * @param other
	 */
	public LatticeBase(LatticeBase other) {
		this.b10 = other.b10;
		this.b11 = other.b11;
		this.b20 = other.b20;
		this.b21 = other.b21;
	}
	
	/**
	 * Reduce this lattice base using the Lagrange-Gauss lattice base reduction algorithm,
	 * following Steven Galbraith "Mathematics of Public Key Cryptography", 2012.
	 * 
	 * @return reduced base
	 */
	public LatticeBase reduce/*_v1*/() {
		// most operations on b1, b2 need long precision, thus we can use long all along
		long b10 = this.b10;
		long b11 = this.b11;
		long b20 = this.b20;
		long b21 = this.b21;
		long tmp0, tmp1;
	
		long B1 = b10*b10 + b11*b11; // ||b1||^2
		
		// b2 = b2 - round(mu)*b1, where mu = inner product(b1, b2) / B1
		int muRounded = (int) (0.5F + (b10*b20 + b11*b21) / (double) B1);
		b20 -= muRounded * b10;
		b21 -= muRounded * b11;
		
		long B2 = b20*b20 + b21*b21; // ||b2||^2
		while (B2 < B1) {
			// swap b1 and b2
			tmp0 = b20;
			tmp1 = b21;
			b20 = b10;
			b21 = b11;
			b10 = tmp0;
			b11 = tmp1;
			
			B1 = B2;
			
			// b2 = b2 - round(mu)*b1, where mu = inner product(b1, b2) / B1
			muRounded = (int) (0.5F + (b10*b20 + b11*b21) / (double) B1);
			b20 -= muRounded * b10;
			b21 -= muRounded * b11;
		
			B2 = b20*b20 + b21*b21; // ||b2||^2
		}
		
		// The result is again in ints; we make Silverman's 'b' positive, which correspond to my b11, b21
		if (b11 < 0) {
			this.b10 = (int) -b10;
			this.b11 = (int) -b11;
		} else {
			this.b10 = (int) b10;
			this.b11 = (int) b11;
		}
		if (b21 < 0) {
			this.b20 = (int) -b20;
			this.b21 = (int) -b21;
		} else {
			this.b20 = (int) b20;
			this.b21 = (int) b21;
		}
		
		// permit chaining
		return this;
	}

	public long determinant() {
		return b10*(long)b21 - b11*(long)b20;
	}
	
	@Override
	public String toString() {
		return "[b10 b20][b11 b21] = [" + b10 + " " + b20 + "][" + b11 + " " + b21 + "]";
	}
}
