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
import java.util.List;

import de.tilman_neumann.jml.factor.base.congruence.AQPair;
import de.tilman_neumann.jml.factor.siqs.tdiv.TDivReport;

public interface TDiv_LQS {

	String getName();

	void initializeForN(int k, BigInteger N, double N_dbl, BigInteger kN, double maxQRest, int[] primesArray, int baseSize, boolean profile);

	void initializeForSpecialQ(int special_q, BQF_xy bqf);

	List<AQPair> testList(ArrayList<IntPair> xyList);

	TDivReport getReport();

	void cleanUp();
}