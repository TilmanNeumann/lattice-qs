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

import java.math.BigInteger;

import org.apache.log4j.Logger;

import de.tilman_neumann.jml.factor.FactorException;
import de.tilman_neumann.jml.modular.JacobiSymbol;
import de.tilman_neumann.jml.primes.exact.AutoExpandingPrimesArray;

/**
 * A SpecialQFinder that sequentially returns small primes having Jacobi(kN, p) == 1), starting with the smallest one.
 * 
 * @author Tilman Neumann
 */
public class SpecialqFinder_small implements SpecialqFinder {
	private static final Logger LOG = Logger.getLogger(SpecialqFinder_small.class);
	private static final boolean DEBUG = false;

	private BigInteger N, kN;
	private int pMax;
	private AutoExpandingPrimesArray primesArray = new AutoExpandingPrimesArray();
	private JacobiSymbol jacobiEngine = new JacobiSymbol();
	private int pIndex;
	
	@Override
	public String getName() {
		return "smallqFinder";
	}

	@Override
	public void initializeForN(BigInteger N, BigInteger kN, int pMax) {
		this.N = N;
		this.kN = kN;
		this.pMax = pMax;
		// compute first pIndex=0
		primesArray.ensureLimit(pMax);
		this.pIndex = 0; // TODO smallest odd?
	}
	
	@Override
	public int nextSpecialQ() throws FactorException {
		while (true) {
			int q = primesArray.getPrime(pIndex);
			pIndex++;
			int jacobi = jacobiEngine.jacobiSymbol(kN, q);
			if (jacobi == 1) {
				if (DEBUG) LOG.info("pMax = " + pMax + " -> specialQ = " + q);
				return q;
			}
			if (jacobi == 0) {
				// Since q is small, this will usually mean that q divides k; but we can try anyway...
				BigInteger qBig = BigInteger.valueOf(q);
				if (N.mod(qBig).intValue()==0) {
					LOG.info("Found factor from special-q = " + q);
					throw new FactorException(qBig);
				}
			}
			// else: jacobi == -1 or 0 but no factor, continue searching
		}		
	}
}
