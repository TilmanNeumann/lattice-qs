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

/**
 * Basic parameters for the quadratic sieve.
 * @author Tilman Neumann
 */
public class LQSSieveParams {
	private static final Logger LOG = Logger.getLogger(LQSSieveParams.class);
	private static final boolean DEBUG = false;

	/** the index of the smallest prime used for sieving. */
	public int pMinIndex;
	/** the smallest prime used for sieving. */
	public int pMin;
	/** the largest prime in the prime base */
	public int pMax;
	/** the area of the full 2-dimensional sieve array */
	public long sieveArrayArea;
	/** maximal Q_rest accepted as smooth candidate */
	public double smoothBound;
	/** sieve array initializer value */
	public byte initializer;
	/** multiplier to scale ln(p) values to the chosen log base */
	public float lnPMultiplier;
	
	public LQSSieveParams(BigInteger kN, int[] primeBase, int primeBaseSize, long sieveArrayArea, int minLnPSumStyle, int specialqSize, double smoothBound, int wantedMinLogPSum) {
		// pMinIndex ~ primeBaseSize^0.33 looks clearly better than primeBaseSize^0.3 or primeBaseSize^0.36.
		// We avoid p[0]==2 which is not used in several sieves.
		this.pMinIndex = Math.max(1, (int) Math.cbrt(primeBaseSize));
		this.pMin = primeBase[pMinIndex];
		this.pMax = primeBase[primeBaseSize-1];
		this.sieveArrayArea = sieveArrayArea;
		this.smoothBound = smoothBound;
		
		// Compute the minimal sum of ln(p_i) values required for a Q to pass the sieve:
		// This bound is computed as average ln(Q/special_q) - ln(smoothBound) values.
		// Assuming a quadratic sieve area, in LQS we have
		// I.   sieve array side length s, sieve array area M = s^2
		// II.  x = e*v1 + f*w1, y = e*v2 + f*w2 with (v1, v2), (w1, w2) being the q-reduced lattice base -> |x|, |y| are of the size of q !!
		// III. Q(x,y) = Ax^2 + 2Bxy + Cy^2 with A ~ sqrt(kN), B ~ sqrt(kN)/2, C ~ -sqrt(kN)
		// Thus we get |Q(x,y)| ~ sqrt(kN) * (1 +- 1 - 1) * (sq)^2 <=> |Q(x,y)| ~ sqrt(kN) * M * q^2
		// => |Q(x,y)|/q ~ sqrt(kN) * M * q. The log of that is ln(kN)/2 + ln(M) + ln(q).
		// Note that this is bigger than |Q(x)| sizes in MPQS/SIQS -> so we should expect LQS to be slower.
		//
		// Interestingly, the wrong case 1 below is fastest and the correct case 3 is the slowest.
		// TODO Is that because case 2, 3 run faster into the cc-problem with 3-partials?
		double minLnPSum;
		switch (minLnPSumStyle) {
		case 1: {
			// Wrong first approach: Subtract ln(special_q)
			// For bigq this requires a smoothBoundExponent ~ 0.1575 // TODO smallq ?
			minLnPSum = Math.log(sieveArrayArea) + Math.log(kN.doubleValue())/2 - Math.log(specialqSize) - Math.log(smoothBound);
			break;
		}
		case 2: {
			// Subtract nothing; just an intermediate choice between the two extremes, case 1 and case 3.
			// For bigq this requires a smoothBoundExponent ~ 0.22
			minLnPSum = Math.log(sieveArrayArea) + Math.log(kN.doubleValue())/2 - Math.log(smoothBound);
			break;
		}
		case 3: {
			// The theoretically correct solution: Add ln(special_q) ~ ln(pMax)
			// For bigq this requires a smoothBoundExponent ~ 0.28
			minLnPSum = Math.log(sieveArrayArea) + Math.log(kN.doubleValue())/2 + Math.log(specialqSize) - Math.log(smoothBound);
			break;
		}
		default: throw new IllegalArgumentException("lnPMinSumStyle must be 1, 2 or 3");
		}
		
		// convert the sieve bound from natural logarithm to the actual logBase:
		float lnLogBase = (float) (minLnPSum / wantedMinLogPSum); // normalizer to be used as a divisor for p_i values
		int minLogPSum = (int) (minLnPSum / lnLogBase); // floor, result should be ~wantedMinLogPSum
		if (DEBUG) {
			float logBase = (float) Math.exp(lnLogBase);
			LOG.debug("logBase=" + logBase + ", lnLogBase=" + lnLogBase + ", minLnPSum = " + minLnPSum + ", minLogPSum = " + minLogPSum);
		}
		initializer = computeInitializerValue(primeBase, pMinIndex, minLogPSum, lnLogBase);
		lnPMultiplier = 1.0F/lnLogBase; // normalizer to be used as a multiplier for p_i values (faster than a divisor)
	}
	
	/**
	 * Compute the initializer value.
	 * @param primesArray prime base
	 * @param pMinIndex the index of the first prime used for sieving
	 * @param minLogPSum
	 * @param lnLogBase
	 * @return initializer byte value
	 */
	private byte computeInitializerValue(int[] primesArray, int pMinIndex, int minLogPSum, double lnLogBase) {
		// compute contribution of small primes in nats
		double lnSmallPSum = 0;
		for (int i=pMinIndex-1; i>=0; i--) {
			int p = primesArray[i];
			lnSmallPSum += Math.log(p) / p;
		}
		// convert value from base e to wanted log base
		double logSmallPSum = lnSmallPSum / lnLogBase;
		// compute initializerValue, rounded
		byte initializerValue = (byte) (128 - minLogPSum + logSmallPSum + 0.5);
		if (DEBUG) LOG.debug("initializerValue = " + initializerValue);
		return initializerValue;
	}

	public byte[] getInitializerBlock() {
		byte[] initializerBlock = new byte[256];
		for (int i=255; i>=0; i--) {
			initializerBlock[i] = initializer;
		}
		return initializerBlock;
	}
}
