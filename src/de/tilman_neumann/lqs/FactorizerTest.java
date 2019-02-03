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

import static de.tilman_neumann.jml.base.BigIntConstants.*;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;

import org.apache.log4j.Logger;

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
	private static final Logger LOG = Logger.getLogger(FactorizerTest.class);

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

			// Trial division: Fastest algorithm for N < 2^29
//			new TDiv31Preload(),
				
			// Lehman: never the best, works until 45 bit
			//new Lehman(),

			// PollardRho:
			// * never the best algorithm
			// * Best BigInteger version is PollardRhoBrent
			//new PollardRho(),
			//new PollardRho_ProductGcd(),
//			new PollardRhoBrent(),
			//new PollardRho31(),

			// SquFoF variants
			// * SquFoF31 is the best algorithm overall for N = 2^29...2^52, SquFoF63 for N = 2^52...2^60
			// * best multiplier sequence = 1680 * {squarefree sequence}
			// * best stopping criterion = O(5.th root(N))
//			new SquFoF63(), // best algorithm for N = 2^52...2^60 (freezes at some N > 2^90)
//			new SquFoF31(), // best algorithm for N = 2^29...2^52
			
			// CFrac
			// * never the best algorithm: SquFoF63 is better for N <= 66 bit, SIQS is better for N >= 55 bits
			// * stopRoot, stopMult: if big enough, then a second k is rarely needed; (5, 1.5) is good
			// * TDiv_CF01 is good for N < 80 bits; for N > 90 bit we need TDiv_CF02
//			new CFrac01(true, 5, 1.5F, 0.152F, 0.253F, new TDiv_CF01(), 10, new MatrixSolver01_Gauss(), 5, false),
//			new CFrac01(true, 5, 1.5F, 0.152F, 0.253F, new TDiv_CF02(), 10, new MatrixSolver01_Gauss(), 5, false),
//			new CFrac63_01(true, 5, 1.5F, 0.152F, 0.25F, new TDiv_CF63_01(), 10, new MatrixSolver01_Gauss(), 12),
//			new CFrac63_01(true, 5, 1.5F, 0.152F, 0.25F, new TDiv_CF63_02(), 10, new MatrixSolver01_Gauss(), 12),

			// SIQS:
			// * TDiv: UBI-variants are better than BigInteger-variants; nLarge is better than 1Large for N >= 200 bit
			// * BlockLanczos is better than Gauss for about N>200 bit
			//new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new SimpleSieve(), new TDiv_QS_1Large(), 10, new MatrixSolver01_Gauss(), false),
//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new Sieve03g(), new TDiv_QS_1Large_UBI(), 10, new MatrixSolver01_Gauss(), false),
//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new Sieve03gU(), new TDiv_QS_1Large_UBI(), 10, new MatrixSolver01_Gauss(), false),

//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new Sieve03g(), new TDiv_QS_1Large_UBI(), 10, new MatrixSolver02_BlockLanczos(), false),
//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new Sieve03g(), new TDiv_QS_2Large_UBI(), 10, new MatrixSolver02_BlockLanczos(), false),
//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new Sieve03g(), new TDiv_QS_nLarge(), 10, new MatrixSolver02_BlockLanczos(), true),
//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new Sieve03g(), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), true),
//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new Sieve03gU(), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), true),

			// sieving with prime powers: best sieve for small N!
//			new SIQS(0.32F, 0.37F, null, null, new PowerOfSmallPrimesFinder(), new SIQSPolyGenerator(), new Sieve03gU(), new TDiv_QS_1Large_UBI(), 10, new MatrixSolver01_Gauss(), false),
//			new SIQS(0.32F, 0.37F, null, null, new PowerOfSmallPrimesFinder(), new SIQSPolyGenerator(), new Sieve03gU(), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), true),
//			new SIQS(0.32F, 0.37F, null, null, new AllPowerFinder(), new SIQSPolyGenerator(), new Sieve03gU(), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), true),
			
			// segmented sieve: slower; best block size is about 32k
//			new SIQS(0.32F, 0.385F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new SingleBlockSieve(32768), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), false),
//			new SIQS(0.32F, 0.385F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new SingleBlockSieveU(32768), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), true),
//			new SIQS(0.32F, 0.41F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new DoubleBlockSieve(32768, 131072), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), false),
	
			// hybrid sieves:
			// * single block hybrid is level with Sieve03g
			// * best block size is about 32k (my processor has 16kB L1 cache, 128kB L2-cache per thread)
//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new SingleBlockHybridSieve01(32768), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), true),
//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new SingleBlockHybridSieveU(32768), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), true),
//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new DoubleBlockHybridSieve01(32768, 131072), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), true),
//			new SIQS(0.32F, 0.37F, null, null, new NoPowerFinder(), new SIQSPolyGenerator(), new DoubleBlockHybridSieveU(32768, 131072), new TDiv_QS_nLarge_UBI(), 10, new MatrixSolver02_BlockLanczos(), true),

			// Multi-threaded SIQS:
			// * 4/6 threads takes over at N around 100 bit (more exact estimates: 4 threads best for N>=88 bits, 6 threads for N>=112 bits)
			// * we need 0.14 < maxQRestExponent < 0.2; everything else is prohibitive; use null for dynamic determination
			// * BlockLanczos is better than Gauss solver for N > 200 bit
//			new PSIQS_U(0.32F, 0.37F, null, null, 6, new NoPowerFinder(), new MatrixSolver02_BlockLanczos(), true),
			new PSIQS_U(0.32F, 0.37F, null, null, 6, new PowerOfSmallPrimesFinder(), new MatrixSolver02_BlockLanczos(), true),
//			new PSIQS_U(0.32F, 0.37F, null, null, 6, new AllPowerFinder(), new MatrixSolver02_BlockLanczos(), true),
//			new PSIQS_SBH_U(0.32F, 0.37F, null, null, 32768, 6, new PowerOfSmallPrimesFinder(), new MatrixSolver02_BlockLanczos(), true), // best for large N

			// Combination of best algorithms for all factor argument sizes
//			new CombinedFactorAlgorithm(4),

			// ====================================================================================================
			// Lattice-QS
			// ====================================================================================================

			// Mmult ~ 0.55 looks best
			// big special q: best so far, has 3-partial problem at "large" N
			/* 3-partial problem at N>=210 bit */ new LQS(0.32F, 0.55F, 1, 0.145F, 10, new SpecialqFinder_big(), new LatticeVectorSieve(), new TDiv_LQS_bigq(), new MatrixSolver02_BlockLanczos(), true),
			///* 3-partial problem for N>=140 bit */ new LQS(0.32F, 0.55F, 2, 0.225F, 10, new SpecialqFinder_big(), new LatticeVectorSieve(), new TDiv_LQS_bigq(), new MatrixSolver02_BlockLanczos(), false),
			///* 3-partial problem for N>=120 bit */ new LQS(0.32F, 0.55F, 3, 0.32F, 10, new SpecialqFinder_big(), new LatticeVectorSieve(), new TDiv_LQS_bigq(), new MatrixSolver02_BlockLanczos(), false),
//			new LQS(0.32F, 0.55F, 1, 0.145F, 10, new SpecialqFinder_big(), new LatticeVectorSieve(), new TDiv_LQS_1L_bigq(), new MatrixSolver02_BlockLanczos(), true),
//			new LQS(0.32F, 0.55F, 2, 0.225F, 10, new SpecialqFinder_big(), new LatticeVectorSieve(), new TDiv_LQS_1L_bigq(), new MatrixSolver02_BlockLanczos(), true),
//			new LQS(0.32F, 0.55F, 3, 0.32F, 10, new SpecialqFinder_big(), new LatticeVectorSieve(), new TDiv_LQS_1L_bigq(), new MatrixSolver02_BlockLanczos(), true),
			// small special q: far too many solver runs
//			new LQS(0.32F, 0.55F, 1, 0.16F, 10, new SpecialqFinder_small(), new LatticeVectorSieve(), new TDiv_LQS_smallq(), new MatrixSolver02_BlockLanczos(), true),
//			new LQS(0.32F, 0.55F, 2, 0.25F, 10, new SpecialqFinder_small(), new LatticeVectorSieve(), new TDiv_LQS_smallq(), new MatrixSolver02_BlockLanczos(), true),
//			new LQS(0.32F, 0.55F, 3, 0.35F, 10, new SpecialqFinder_small(), new LatticeVectorSieve(), new TDiv_LQS_smallq(), new MatrixSolver02_BlockLanczos(), true),
//			new LQS_smallq(0.32F, 0.55F, 1, 0.22F, 10, new SpecialqFinder_small(), new LatticeVectorSieve(), new TDiv_LQS_smallq(), new MatrixSolver02_BlockLanczos(), true),
//			new LQS_smallq(0.32F, 0.55F, 2, 0.26F, 10, new SpecialqFinder_small(), new LatticeVectorSieve(), new TDiv_LQS_smallq(), new MatrixSolver02_BlockLanczos(), true),
//			new LQS_smallq(0.32F, 0.55F, 3, 0.285F, 10, new SpecialqFinder_small(), new LatticeVectorSieve(), new TDiv_LQS_smallq(), new MatrixSolver02_BlockLanczos(), true)
		};
	}
	
	private void testRange(int bits) {
		BigInteger N_min = I_1.shiftLeft(bits-1);
		// find N-set for square tests
		//ArrayList NSet = TestsetGenerator.generate(bits, N_COUNT);
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
