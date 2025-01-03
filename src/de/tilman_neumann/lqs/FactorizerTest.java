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

import static de.tilman_neumann.jml.base.BigIntConstants.*;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;

import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

import de.tilman_neumann.jml.factor.CombinedFactorAlgorithm;
import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.factor.TestNumberNature;
import de.tilman_neumann.jml.factor.TestsetGenerator;
import de.tilman_neumann.jml.factor.base.congruence.*;
import de.tilman_neumann.jml.factor.base.matrixSolver.*;
import de.tilman_neumann.jml.factor.cfrac.*;
import de.tilman_neumann.jml.factor.cfrac.tdiv.*;
import de.tilman_neumann.jml.factor.lehman.*;
import de.tilman_neumann.jml.factor.pollardRho.*;
import de.tilman_neumann.jml.factor.psiqs.*;
import de.tilman_neumann.jml.factor.siqs.*;
import de.tilman_neumann.jml.factor.siqs.poly.SIQSPolyGenerator;
import de.tilman_neumann.jml.factor.siqs.poly.baseFilter.*;
import de.tilman_neumann.jml.factor.siqs.powers.*;
import de.tilman_neumann.jml.factor.siqs.sieve.*;
import de.tilman_neumann.jml.factor.siqs.tdiv.*;
import de.tilman_neumann.jml.factor.squfof.*;
import de.tilman_neumann.jml.factor.tdiv.*;
import de.tilman_neumann.jml.sequence.*;
import de.tilman_neumann.util.ConfigUtil;
import de.tilman_neumann.util.TimeUtil;

/**
 * Main class to compare the performance of factor algorithms; here, in particular with the lattice QS algorithm.
 * @author Tilman Neumann
 */
@SuppressWarnings("unused") // suppress warnings on unused imports
public class FactorizerTest {
	private static final Logger LOG = LogManager.getLogger(FactorizerTest.class);

	// algorithm options
	/** number of test numbers */
	private static final int N_COUNT = 1;
	/** the bit size of N to start with */
	private static final int START_BITS = 70;
	/** the increment in bit size from test set to test set */
	private static final int INCR_BITS = 10;
	/** maximum number of bits to test (no maximum if null) */
	private static final Integer MAX_BITS = null;
	/** each algorithm is run REPEATS times for each input in order to reduce GC influence on timings */
	private static final int REPEATS = 1;

	/** 
	 * Algorithms to compare. Non-static to permit to use Loggers in the algorithm constructors.
	 */
	private FactorAlgorithm[] algorithms;
	
	public FactorizerTest() {
		algorithms = new FactorAlgorithm[] {

			// CFrac
//			new CFrac(true, 5, 1.5F, 0.152F, 0.253F, new TDiv_CF02(), 10, new MatrixSolverGauss02(), 5),

			// SIQS
			new SIQS(0.32F, 0.37F, null, new NoPowerFinder(), new SIQSPolyGenerator(), new Sieve03gU(), new TDiv_QS_nLP_Full(true), 10, new MatrixSolverBlockLanczos()),

			// Multi-threaded SIQS
//			new PSIQS_U(0.32F, 0.37F, null, 6, new NoPowerFinder(), new MatrixSolverBlockLanczos()),

			// Combination of best algorithms for all factor argument sizes
//			new CombinedFactorAlgorithm(4),

			// ====================================================================================================
			// Lattice-QS
			// ====================================================================================================

			// Mmult ~ 0.55 looks best
			// big special q: best so far, has 3-partial problem at "large" N
			/* 3-partial problem at N>=210 bit */ new LQS(0.32F, 0.55F, 1, 0.145F, 10, new SpecialqFinder_big(), new LatticeVectorSieve(), new TDiv_LQS_bigq(), new MatrixSolverBlockLanczos()),
//			/* 3-partial problem for N>=140 bit */ new LQS(0.32F, 0.55F, 2, 0.225F, 10, new SpecialqFinder_big(), new LatticeVectorSieve(), new TDiv_LQS_bigq(), new MatrixSolverBlockLanczos()),
//			/* 3-partial problem for N>=120 bit */ new LQS(0.32F, 0.55F, 3, 0.32F, 10, new SpecialqFinder_big(), new LatticeVectorSieve(), new TDiv_LQS_bigq(), new MatrixSolverBlockLanczos()),
//			new LQS(0.32F, 0.55F, 1, 0.145F, 10, new SpecialqFinder_big(), new LatticeVectorSieve(), new TDiv_LQS_1L_bigq(), new MatrixSolverBlockLanczos()),
//			new LQS(0.32F, 0.55F, 2, 0.225F, 10, new SpecialqFinder_big(), new LatticeVectorSieve(), new TDiv_LQS_1L_bigq(), new MatrixSolverBlockLanczos()),
//			new LQS(0.32F, 0.55F, 3, 0.32F, 10, new SpecialqFinder_big(), new LatticeVectorSieve(), new TDiv_LQS_1L_bigq(), new MatrixSolverBlockLanczos()),
			// small special q: far too many solver runs
//			new LQS(0.32F, 0.55F, 1, 0.16F, 10, new SpecialqFinder_small(), new LatticeVectorSieve(), new TDiv_LQS_smallq(), new MatrixSolverBlockLanczos()),
//			new LQS(0.32F, 0.55F, 2, 0.25F, 10, new SpecialqFinder_small(), new LatticeVectorSieve(), new TDiv_LQS_smallq(), new MatrixSolverBlockLanczos()),
//			new LQS(0.32F, 0.55F, 3, 0.35F, 10, new SpecialqFinder_small(), new LatticeVectorSieve(), new TDiv_LQS_smallq(), new MatrixSolverBlockLanczos()),
//			new LQS_smallq(0.32F, 0.55F, 1, 0.22F, 10, new SpecialqFinder_small(), new LatticeVectorSieve(), new TDiv_LQS_smallq(), new MatrixSolverBlockLanczos()),
//			new LQS_smallq(0.32F, 0.55F, 2, 0.26F, 10, new SpecialqFinder_small(), new LatticeVectorSieve(), new TDiv_LQS_smallq(), new MatrixSolverBlockLanczos()),
//			new LQS_smallq(0.32F, 0.55F, 3, 0.285F, 10, new SpecialqFinder_small(), new LatticeVectorSieve(), new TDiv_LQS_smallq(), new MatrixSolverBlockLanczos())
		};
	}
	
	private void testRange(int bits) {
		BigInteger N_min = I_1.shiftLeft(bits-1);
		BigInteger[] testNumbers = TestsetGenerator.generate(N_COUNT, bits, TestNumberNature.MODERATE_SEMIPRIMES);
		LOG.info("Test N with " + bits + " bits, i.e. N >= " + N_min);
		
		// take 3 timings for each algorithm to be quite sure that one timing is not falsified by garbage collection
		TreeMap<Long, List<FactorAlgorithm>> ms_2_algorithms = new TreeMap<Long, List<FactorAlgorithm>>();
		for (int i=0; i<REPEATS; i++) {
			for (FactorAlgorithm algorithm : algorithms) {
				// exclude special size implementations
				String algName = algorithm.getName();
				if (bits<45 && algName.startsWith("SIQS")) continue; // unstable for smaller N
				if (bits<57 && algName.startsWith("PSIQS")) continue; // unstable for smaller N
				if (bits>98 && algName.startsWith("CFrac63")) continue; // unstable for N>98 bits
				if (bits>63 && algName.startsWith("TDiv63")) continue; // long implementation
				if (bits>52 && algName.equals("SquFoF31")) continue; // int implementation
				if (bits>45 && algName.startsWith("Lehman")) continue; // int implementation
				if (bits>31 && algName.startsWith("TDiv31")) continue; // int implementation
				if (bits>31 && algName.startsWith("PollardRho31")) continue; // long implementation
				
				System.gc(); // create equal conditions for all algorithms

				int failCount = 0;
				long startTimeMillis = System.currentTimeMillis();
				for (BigInteger N : testNumbers) {
					BigInteger factor = algorithm.findSingleFactor(N);
					// test correctness
					if (factor==null || factor.equals(I_0) || factor.equals(I_1) || factor.mod(N).equals(I_0)) {
						//LOG.error("FactorAlgorithm " + algorithm.getName() + " did not find a factor of N=" + N + ", it returned " + factor);
						failCount++;
					} else {
						// not null, not trivial -> test division
						BigInteger[] test = N.divideAndRemainder(factor);
						if (!test[1].equals(I_0)) {
							//LOG.error("FactorAlgorithm " + algorithm.getName() + " returned " + factor + ", but this is not a factor of N=" + N);
							failCount++;
						}
					}
				}
				long endTimeMillis = System.currentTimeMillis();
				long duration = endTimeMillis - startTimeMillis; // duration in ms
				//LOG.debug("algorithm " + algName + " finished test set with " + bits + " bits");
				List<FactorAlgorithm> algList = ms_2_algorithms.get(duration);
				if (algList==null) algList = new ArrayList<FactorAlgorithm>();
				algList.add(algorithm);
				ms_2_algorithms.put(duration, algList);
				if (failCount>0) {
					LOG.error("FactorAlgorithm " + algorithm.getName() + " failed at " + failCount + "/" + N_COUNT + " test numbers...");
				}
			}
		}
		
		// log best algorithms first
		int rank=1;
		for (long ms : ms_2_algorithms.keySet()) {
			List<FactorAlgorithm> algList = ms_2_algorithms.get(ms);
			int j=0;
			for (FactorAlgorithm algorithm : algList) {
				String durationStr = TimeUtil.timeStr(ms);
				LOG.info("#" + rank + ": Algorithm " + algorithm.getName() + " took " + durationStr);
				j++;
			}
			rank += j;
		}
	}

	/**
	 * Test factor algorithms for sets of factor arguments of growing size and report timings after each set.
	 * @param args ignored
	 */
	public static void main(String[] args) {
    	ConfigUtil.initProject();
    	FactorizerTest testEngine = new FactorizerTest();
		int bits = START_BITS;
		while (true) {
			// test N with the given number of bits, i.e. 2^(bits-1) <= N <= (2^bits)-1
			testEngine.testRange(bits);
			bits += INCR_BITS;
			if (MAX_BITS!=null && bits > MAX_BITS) break;
		}
	}
}
