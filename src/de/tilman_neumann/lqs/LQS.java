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

import static de.tilman_neumann.jml.factor.base.GlobalFactoringOptions.*;
import static de.tilman_neumann.jml.base.BigIntConstants.*;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

import de.tilman_neumann.jml.factor.FactorAlgorithm;
import de.tilman_neumann.jml.factor.FactorException;
import de.tilman_neumann.jml.factor.base.PrimeBaseGenerator;
import de.tilman_neumann.jml.factor.base.congruence.AQPair;
import de.tilman_neumann.jml.factor.base.congruence.CongruenceCollector;
import de.tilman_neumann.jml.factor.base.congruence.CongruenceCollector01;
import de.tilman_neumann.jml.factor.base.congruence.CongruenceCollectorReport;
import de.tilman_neumann.jml.factor.base.matrixSolver.FactorTest;
import de.tilman_neumann.jml.factor.base.matrixSolver.FactorTest01;
import de.tilman_neumann.jml.factor.base.matrixSolver.MatrixSolver;
import de.tilman_neumann.jml.factor.base.matrixSolver.MatrixSolver_BlockLanczos;
import de.tilman_neumann.jml.factor.siqs.KnuthSchroeppel;
import de.tilman_neumann.jml.factor.siqs.ModularSqrtsEngine;
import de.tilman_neumann.jml.factor.siqs.sieve.SieveReport;
import de.tilman_neumann.jml.factor.siqs.tdiv.TDivReport;
import de.tilman_neumann.jml.powers.PurePowerTest;
import de.tilman_neumann.jml.primes.exact.AutoExpandingPrimesArray;
import de.tilman_neumann.util.ConfigUtil;
import de.tilman_neumann.util.TimeUtil;
import de.tilman_neumann.util.Timer;

/**
 * Lattice QS, following an idea by R.S.Silverman in http://www.mersenneforum.org/showthread.php?t=14080.
 * 
 * So far, the biggest number factored with this algorithm had 220 bits.
 * 
 * @author Tilman Neumann
 */
public class LQS extends FactorAlgorithm {
	private static final Logger LOG = LogManager.getLogger(LQS.class);
	private static final boolean DEBUG = false;

	private PurePowerTest powerTest = new PurePowerTest();
	private KnuthSchroeppel multiplierFinder = new KnuthSchroeppel(); // used to compute the multiplier k
	// TODO use KS without correction for d=2
	// TODO do I have to check 4kN for a rich prime base instead of kN? Or is it the same because 4=2^2 is square?
	
	// prime base
	private float Cmult; // multiplier to compute prime base size
	private PrimeBaseGenerator primeBaseBuilder = new PrimeBaseGenerator();
	private AutoExpandingPrimesArray primesArray = AutoExpandingPrimesArray.get();
	
	private ModularSqrtsEngine modularSqrtsEngine = new ModularSqrtsEngine(); // computes tArray

	// poly
	private SpecialqFinder specialqFinder;
	
	// sieve
	private LatticeSieve latticeSieve;
	private float Mmult;
	private int minLnPSumStyle;
	private Float smoothBoundExponent0;
	private float smoothBoundExponent;
	private int sieveArraySideLength;
	
	// trial division engine
	private TDiv_LQS auxFactorizer;
	// extra congruences to have a bigger chance that the equation system solves. the likelihood is >= 1-2^(extraCongruences+1)
	private int extraCongruences;
	
	// collects the congruences we find
	private CongruenceCollector congruenceCollector;
	// The solver used for smooth congruence equation systems.
	private MatrixSolver matrixSolver;
	
	// statistics
	private Timer timer = new Timer();
	private long powerTestDuration, initNDuration, initPolyDuration;

	public LQS(float Cmult, float Mmult, int minLnPSumStyle, Float smoothBoundExponent, int extraCongruences, SpecialqFinder specialqFinder, LatticeSieve latticeSieve, TDiv_LQS tdivEngine, MatrixSolver matrixSolver) {
		this.Cmult = Cmult;
		this.Mmult = Mmult;
		this.minLnPSumStyle = minLnPSumStyle;
		this.smoothBoundExponent0 = smoothBoundExponent;
		this.specialqFinder = specialqFinder;
		this.latticeSieve = latticeSieve;
		this.auxFactorizer = tdivEngine;
		this.extraCongruences = extraCongruences;
		this.congruenceCollector = new CongruenceCollector01();
		this.matrixSolver = matrixSolver;
	}

	@Override
	public String getName() {
		String smoothBoundExponentStr = "smoothBound=" + String.format("%.3f", smoothBoundExponent);
		return "LQS(Cmult=" + Cmult + ", Mmult=" + Mmult + ", sieveArraySideLength=" + sieveArraySideLength + ", minLnPSumStyle = " + minLnPSumStyle + ", " + smoothBoundExponentStr + ", " + specialqFinder.getName() + ", " + latticeSieve.getName() + ", " + auxFactorizer.getName() + ", " + matrixSolver.getName() + ")";
	}
	
	/**
	 * Test the current N.
	 * @return factor, or null if no factor was found.
	 */
	public BigInteger findSingleFactor(BigInteger N) {
		if (ANALYZE) {
			timer.start();
			powerTestDuration = initNDuration = initPolyDuration = 0;
		}

		// the quadratic sieve does not work for pure powers; check that first:
		PurePowerTest.Result purePower = powerTest.test(N);
		if (purePower!=null) {
			// N is indeed a pure power -> return a factor that is about sqrt(N)
			BigInteger factor = purePower.base.pow(purePower.exponent>>1);
			LOG.info("N is a pure power -> found factor " + factor);
			return factor;
		} // else: no pure power, run quadratic sieve
		if (ANALYZE) powerTestDuration += timer.capture();

		// compute Knuth-Schroppel multiplier
		int k = multiplierFinder.computeMultiplier(N);
		BigInteger kN = N; // iff k==1
		if (k>1) {
			BigInteger kBig = BigInteger.valueOf(k);
			// avoid square kN without square N; that would lead to an infinite loop in trial division
			if (N.mod(kBig).equals(I_0)) {
				LOG.info("k = " + k + " divides N -> Found factor.");
				return kBig;
			}
			kN = kBig.multiply(N);
		}
		
		// Compute prime base size:
		// http://www.mersenneforum.org/showthread.php?s=3087aa210d8d7f1852c690a45f22d2e5&t=11116&page=2:
		// "A rough estimate of the largest prime in the factor base is exp(0.5 * sqrt(ln(N) * ln(ln(N))))".
		// Here we estimate the number of entries instead of the largest element in the prime base,
		// with Cmult instead of the constant 0.5. A good value is Cmult = 0.32.
		int NBits = N.bitLength();
		double N_dbl = N.doubleValue();
		double lnN = Math.log(N_dbl);
		double lnTerm = Math.sqrt(lnN * Math.log(lnN)); // (lnN)^0.5 * (lnlnN)^(1-0.5)
		double primeBaseSize_dbl = Math.exp(Cmult * lnTerm);
		if (primeBaseSize_dbl > Integer.MAX_VALUE) {
			LOG.error("N=" + N + " (" + NBits + " bits) is too big for LQS!");
			return null;
		}
		int primeBaseSize = Math.max(30, (int) primeBaseSize_dbl); // min. size for very small N

		// The number of congruences we need to find before we try to solve the smooth congruence equation system:
		// We want: #equations = #variables + some extra congruences
		int requiredSmoothCongruenceCount = primeBaseSize + extraCongruences;
		if (DEBUG) LOG.info("required number of smooths = " + requiredSmoothCongruenceCount);
		
		// Create the reduced prime base for kN
		int[] primesArray = new int[primeBaseSize];
		primeBaseBuilder.computeReducedPrimeBase(kN, primeBaseSize, primesArray);
		int pMax = primesArray[primeBaseSize-1];
		if (DEBUG) LOG.debug("Prime base = " + Arrays.toString(primesArray));

		// Compute the t with t^2 == kN (mod p) for all p: Throws a FactorException if some p divides N
		int[] tArray = modularSqrtsEngine.computeTArray(primesArray, primeBaseSize, kN);

		// Compute the area of the 2-dim. sieve array:
		// We assume it is similar to the 1D sieve array size in ordinary QS.
		// We want the side length of the sieve array to be even to have integer half side length.
		long sieveArrayArea = (int) (6144 + Math.exp(Mmult * lnTerm)); // 6144 = best experimental result for small N
		sieveArraySideLength = (int) Math.sqrt(sieveArrayArea);
		if ((sieveArraySideLength & 1)==1) sieveArraySideLength+=1;
		if (DEBUG) LOG.debug("N=" + N + ", k=" + k + ": sieveArrayArea = " + sieveArrayArea + ", sieveArraySideLength = " + sieveArraySideLength);
		
		// Construct a binary quadratic form Q(x, y) = A*x^2 + 2Bxy + C*y^2 with discriminant D=4kN.
		BQF_xy bqf = null;
		try {
			bqf = new BQFFactory().createBQFWithDiscriminant4kN(k, N, kN, primesArray, primeBaseSize, sieveArrayArea);
			if (DEBUG) LOG.debug(bqf.toString());
		} catch (FactorException fe) {
			BigInteger factor = fe.getFactor();
			LOG.info("Found factor " + factor + " of " + kN + " during BQF construction!");
			return fe.getFactor();
		}
		
		// compute biggest QRest admitted for a smooth relation
		if (smoothBoundExponent0 != null) {
			smoothBoundExponent = smoothBoundExponent0;
		} else {
			smoothBoundExponent = (NBits<=150) ? 0.16F : 0.16F + (NBits-150.0F)/5250;
		}
		double smoothBound = Math.pow(N_dbl, smoothBoundExponent);

		// compute some basic parameters for N
		int specialqSize = pMax;
		LQSSieveParams sieveParams = new LQSSieveParams(kN, primesArray, primeBaseSize, sieveArrayArea, minLnPSumStyle, specialqSize, smoothBound, 127);
		// compute logP array
		byte[] logPArray = computeLogPArray(primesArray, primeBaseSize, sieveParams.lnPMultiplier);
		
		// Allocate 2-dimensional sieve array
		byte[][] sieveArray = new byte[sieveArraySideLength][sieveArraySideLength];
		byte[] initializedSieveLine = new byte[sieveArraySideLength];
		byte initializer = sieveParams.initializer;
		for (int i=sieveArraySideLength-1; i>=0; i--) {
			initializedSieveLine[i] = initializer;
		}
		
		// dontUseArray[e][f]==1 means that we don't want to use (e,f) for sieving.
		byte[][] dontUseArray = computeDontUseArray(sieveArraySideLength);
		
		// initialize sub-algorithms for new N
		latticeSieve.initializeForN(k, N, kN, primesArray, tArray, primeBaseSize, sieveParams, bqf);
		auxFactorizer.initializeForN(k, N, N_dbl, kN, smoothBound, primesArray, primeBaseSize);
		FactorTest factorTest = new FactorTest01(N);
		matrixSolver.initialize(N, factorTest);
		congruenceCollector.initialize(N, primeBaseSize, matrixSolver, factorTest);

		specialqFinder.initializeForN(N, kN, pMax);
		if (ANALYZE) initNDuration += timer.capture();
			
		// "special_q" loop
		BigInteger factor = null;
		while (factor == null) {
			factor = testSpecialQ(bqf, logPArray, sieveArray, initializedSieveLine, dontUseArray);
		}
		
		if (ANALYZE) {
			// get all reports
			SieveReport sieveReport = latticeSieve.getReport();
			TDivReport tdivReport = auxFactorizer.getReport();
			CongruenceCollectorReport ccReport = congruenceCollector.getReport();
			// solverReport is not urgently needed
			
			long sieveDuration = sieveReport.getTotalDuration(1);
			long tdivDuration = tdivReport.getTotalDuration(1);
			
			// report results
			LOG.info(getName() + ":");
			LOG.info("Found factor " + factor + " (" + factor.bitLength() + " bits) of N=" + N + " (" + NBits + " bits) in " + TimeUtil.timeStr(timer.totalRuntime()));
			int pMaxBits = 32 - Integer.numberOfLeadingZeros(pMax);
			LOG.info("    multiplier k = " + k + ", kN%8 = " + kN.mod(I_8) + ", primeBaseSize = " + primeBaseSize + ", pMax = " + pMax + " (" + pMaxBits + " bits), sieveArrayArea = " + sieveArrayArea);
			LOG.info("    tDiv: " + tdivReport.getOperationDetails());
			LOG.info("    cc: " + ccReport.getOperationDetails());
			if (ANALYZE_LARGE_FACTOR_SIZES) {
				for (int i=1; i<=2; i++) LOG.info("        " + ccReport.getSmoothBigFactorPercentiles(i));
				for (int i=1; i<=2; i++) LOG.info("        " + ccReport.getSmoothQRestPercentiles(i));
				for (int i=1; i<=2; i++) LOG.info("        " + ccReport.getPartialBigFactorPercentiles(i));
				for (int i=1; i<=2; i++) LOG.info("        " + ccReport.getPartialQRestPercentiles(i));
				LOG.info("        " + ccReport.getNonIntFactorPercentages());
			}
			if (ANALYZE_Q_SIGNS) {
				LOG.info("        " + ccReport.getPartialQSignCounts());
				LOG.info("        " + ccReport.getSmoothQSignCounts());
			}
			LOG.info("    #solverRuns = " + congruenceCollector.getSolverRunCount() + ", #tested null vectors = " + congruenceCollector.getTestedNullVectorCount());
			LOG.info("    Approximate phase timings: powerTest=" + powerTestDuration + "ms, initN=" + initNDuration + "ms, initPoly=" + initPolyDuration + "ms, sieve=" + sieveDuration + "ms, tdiv=" + tdivDuration + "ms, cc=" + congruenceCollector.getCollectDuration() + "ms, solver=" + congruenceCollector.getSolverDuration() + "ms");
			LOG.info("    -> sieve sub-timings: " + sieveReport.getPhaseTimings(1));
			
			// TDiv, CC and solver have no sub-timings yet
		}
		// TODO clean up
		
		// return factor
		return factor;
	}
	
	private byte[] computeLogPArray(int[] primesArray, int primeBaseSize, float lnPMultiplier) {
		byte[] logPArray = new byte[primeBaseSize];
		for (int i=primeBaseSize-1; i>=0; i--) {
			logPArray[i] = (byte) ((float) Math.log(primesArray[i]) * lnPMultiplier + 0.5F);
		}
		return logPArray;
	}
	
	private byte[][] computeDontUseArray(int sieveArraySideLength) {
		// dontUseArray[e][f]==1 means that we don't want to use (e,f) for sieving.
		// We do not want (e,f)=(0,0) and those with gcd(e,f)>1.
		// We must take into account that e runs from -halfSieveArraySideLength to halfSieveArraySideLength-1.
		byte[][] dontUseArray = new byte[sieveArraySideLength][sieveArraySideLength];
		int halfSieveArraySideLength = sieveArraySideLength>>1;
		dontUseArray[halfSieveArraySideLength][0] = 1;
		for (int i=0; ; i++) {
			int p = primesArray.getPrime(i);
			if (p>=sieveArraySideLength) break;
			if (DEBUG) LOG.debug("p_" + i + " = " + p);
			for (int e=p; e<halfSieveArraySideLength; e+=p) {
				for (int f=p; f<sieveArraySideLength; f+=p) {
					dontUseArray[e+halfSieveArraySideLength][f] = 1;
				}
			}
			for (int e=-p; e>=-halfSieveArraySideLength; e-=p) {
				for (int f=p; f<sieveArraySideLength; f+=p) {
					dontUseArray[e+halfSieveArraySideLength][f] = 1;
				}
			}
		}
		if (DEBUG) {
			LOG.debug("dontUseArray = ");
			for (int f=0; f<sieveArraySideLength; f++) {
				String rowStr = "f=" + f + ": ";
				for (int e=0; e<sieveArraySideLength; e++) {
					rowStr += dontUseArray[e][f];
				}
				LOG.debug(rowStr);
			}
		}
		return dontUseArray;
	}
	
	private BigInteger testSpecialQ(BQF_xy bqf, byte[] logPArray, byte[][] sieveArray, byte[] initializedSieveLine, byte[][] dontUseArray) {
		if (ANALYZE) timer.capture();
		// Choose some "special_q" generating the lattice L of Q(x, y) divisible by q.
		// We want q to be int, and slightly bigger than the prime base.
		int specialQ;
		try {
			specialQ = specialqFinder.nextSpecialQ();
		} catch (FactorException fe) {
			return fe.getFactor();
		}
		// TODO Perform tests to find the optimal size of special q. Also consider q < pMax !
	
		auxFactorizer.initializeForSpecialQ(specialQ, bqf);
		if (ANALYZE) initPolyDuration += timer.capture();

		// do the lattice sieve
		ArrayList<IntPair> smoothCandidatesFromLatticeSieve = latticeSieve.sieve(specialQ, logPArray, sieveArray, initializedSieveLine, dontUseArray, sieveArraySideLength);
		if (DEBUG) if (smoothCandidatesFromLatticeSieve.size() > 0) LOG.debug("q=" + specialQ + ": Lattice sieve found " + smoothCandidatesFromLatticeSieve.size() + " smooth candidates");

		// trial division stage: produce AQ-pairs
		List<AQPair> aqPairs = this.auxFactorizer.testList(smoothCandidatesFromLatticeSieve);
		if (DEBUG) LOG.debug("Trial division found " + aqPairs.size() + " Q(x) smooth enough for a congruence.");

		// add all congruences and run matrix solver if appropriate
		congruenceCollector.collectAndProcessAQPairs(aqPairs);
		return congruenceCollector.getFactor();
	}

	// Standalone test --------------------------------------------------------------------------------------------------

	/**
	 * Easy test numbers
	 * 11111111111111111111111111
	 * 5679148659138759837165981543
	 * 1111111111111111111111111111711 = 5689 * 172657 * 549540617 * 2058436950671
	 * 11111111111111111111111111117777711
	 * 111111111111111111111111111177777111
	 * 11111111111111111111111111117777711111
	 * 111111111111111111111111111177777111111
	 * 1111111111111111111111111111777771111111
	 * 11111111111111111111111111117777711111113
	 * 11111111111111111111111111155555555555111111111111111
	 * 
	 * (small but) hard semiprimes:
	 * 240064589223751 = 15485863 * 15502177
	 * 240316575186487 = 15485863 * 15518449
	 * 240569743163473 = 15502177 * 15518449
	 * 646812682155893468167 = 1858312369 * 348064562743
	 * 1524157998833865801420881 = 1234567890133 * 1234567990157
	 * 1020015712220023286172890929 = 2064262302587 * 494130862604867
	 * 
	 * "Too many solver runs" fixed:
	 * 689396972318681563762687 = 96152311 * 7169842982959417
	 * 408013346975910849946449119300452669 = 1213888521319 * 336120936815983179646651
	 * 880295022091761685646642583312545957282042227 = 918932272595419 * 957954191341511139128076418633
	 * 
	 * TODO Too many solver runs:
	 * 733290673195330967156941 = 1782090407 * 411477818585963 (4 runs)
	 * 
	 * TODO Infinite loop ?
	 * 11111111111111111111111111117 prime
	 * 11111111111111111111111111117777711111111 prime
	 * 1114170342238328228419686817856626026481517749572426317489981 (running 50 minutes instead of expected 8-12)
	 */
	private static void testInput() {
		LQS qs = new LQS(0.32F, 0.55F, 1, null, 10, new SpecialqFinder_big(), new LatticeVectorSieve(), new TDiv_LQS_bigq(), new MatrixSolver_BlockLanczos());
		//LQS qs = new LQS(0.32F, 0.55F, 1, null, 10, new SpecialqFinder_small(), new LatticeVectorSieve(), new TDiv_LQS_smallq(), new MatrixSolver02_BlockLanczos(), true);

		Timer timer = new Timer();
		while(true) {
			try {
				LOG.info("Please insert the number to factor:");
				BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
				String line = in.readLine();
				String input = line != null ? line.trim() : "";
				//LOG.debug("input = >" + input + "<");
				BigInteger N = new BigInteger(input);
				LOG.info("Factoring " + N + " (" + N.bitLength() + " bits)...");
				timer.capture();
				BigInteger factor = qs.findSingleFactor(N);
				if (factor != null) {
					long duration = timer.capture();
					LOG.info("Found factor " + factor + " in " + TimeUtil.timeStr(duration) + ".");
				} else {
					LOG.info("No factor found...");
				}
			} catch (Exception ex) {
				LOG.error("Error: " + ex, ex);
			}
		}
	}
	
	/**
	 * Test of input k, N and #iterations.
	 * @param args ignored
	 */
	public static void main(String[] args) {
    	ConfigUtil.initProject();
    	testInput();
	}
}
