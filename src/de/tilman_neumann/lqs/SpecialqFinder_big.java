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

import de.tilman_neumann.jml.factor.FactorException;
import de.tilman_neumann.jml.modular.JacobiSymbol;
import de.tilman_neumann.jml.primes.exact.AutoExpandingPrimesArray;

/**
 * A SpecialQFinder that sequentially returns the primes > pMax having Jacobi(kN, p) == 1.
 * 
 * @author Tilman Neumann
 */
public class SpecialqFinder_big implements SpecialqFinder {
	private static final Logger LOG = Logger.getLogger(SpecialqFinder_big.class);
	private static final boolean DEBUG = false;

	private BigInteger N, kN;
	private int pMax;
	private AutoExpandingPrimesArray primesArray = new AutoExpandingPrimesArray();
	private JacobiSymbol jacobiEngine = new JacobiSymbol();
	private int pIndex;
	
	@Override
	public String getName() {
		return "bigqFinder";
	}
	
	@Override
	public void initializeForN(BigInteger N, BigInteger kN, int pMax) {
		this.N = N;
		this.kN = kN;
		this.pMax = pMax;
		// compute first pIndex as the index of the first prime > pMax
		this.pIndex = primesArray.ensureLimit(2*pMax).getInsertPosition(pMax);
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
				// Since usually k << q, this is likely a factor of N!
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
