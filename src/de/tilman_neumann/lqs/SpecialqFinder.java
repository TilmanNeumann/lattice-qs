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

import de.tilman_neumann.jml.factor.FactorException;

/**
 * Interface for special_q selection algorithms.
 * Special q must satisfy Jacobi(kN, q) > 1; otherwise polynomial roots computation in class BQF would be slightly wrong.
 * 
 * @author Tilman Neumann
 */
public interface SpecialqFinder {
	void initializeForN(BigInteger N, BigInteger kN, int pMax);
	
	String getName();
	
	/**
	 * @return special_q
	 * @throws FactorException 
	 */
	int nextSpecialQ() throws FactorException;
}
