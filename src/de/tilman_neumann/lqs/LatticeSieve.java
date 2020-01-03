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
import java.util.ArrayList;

import de.tilman_neumann.jml.factor.FactorException;
import de.tilman_neumann.jml.factor.siqs.sieve.SieveReport;

/**
 * Interface for 2D lattice sievers.
 * @author Tilman Neumann
 */
public interface LatticeSieve {

	String getName();
	
	void initializeForN(int k, BigInteger N, BigInteger kN, int[] primesArray, int[] tArray, int primeBaseSize, LQSSieveParams sieveParams, BQF_xy Qxy);
	
	ArrayList<IntPair> sieve(int q, byte[] logPArray, byte[][] sieveArray, byte[] initializedSieveLine, byte[][] dontUseArray, int sieveArraySideLength) throws FactorException;

	public SieveReport getReport();

}
