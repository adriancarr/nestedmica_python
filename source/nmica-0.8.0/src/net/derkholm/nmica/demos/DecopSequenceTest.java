/*
 * NestedMICA Motif Inference Toolkit
 *
 * Copyright (c) 2004-2007: Genome Research Ltd.
 *
 * NestedMICA is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * or see the on-line version at http://www.gnu.org/copyleft/lgpl.txt
 *
 */
/*
 * Created on Sep 19, 2004
 */
package net.derkholm.nmica.demos;

import net.derkholm.nmica.utils.ArrayTools;
import net.derkholm.nmica.utils.CollectTools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringTokenizer;

/**
 * @author thomas
 */
public class DecopSequenceTest {
    public static void main(String[] args) 
    		throws Exception
    {
        File scoreListFile = new File(args[0]);
        double target = Double.parseDouble(args[1]);
        
        List<Double> negScoreList = new ArrayList<Double>();
        List<Double> posScoreList = new ArrayList<Double>();
        
        BufferedReader br = new BufferedReader(new FileReader(scoreListFile));
        for (String line = br.readLine(); line != null; line = br.readLine()) {
            StringTokenizer toke = new StringTokenizer(line);
            negScoreList.add(new Double(toke.nextToken()));
            posScoreList.add(new Double(toke.nextToken()));
        }
        
        double[] negScores = CollectTools.toDoubleArray(negScoreList);
        double[] posScores = CollectTools.toDoubleArray(posScoreList);
        double[] difScores = new double[negScores.length];
        
        double negScore = 0;
        double posScore = 0;
        double optScore = 0;
        for (int d = 0; d < negScores.length; ++d) {
            difScores[d] = posScores[d] - negScores[d];
            negScore += negScores[d];
            posScore += posScores[d];
            optScore += Math.max(negScores[d], posScores[d]);
        }
        
        System.err.println("Everything off: " + negScore);
        System.err.println("Everything on: " + posScore);
        System.err.println("Optimum: " + optScore);
        
        double[] difScoresSorted = (double[]) ArrayTools.copy(difScores);
        Arrays.sort(difScoresSorted);
        
        double currentScore = negScore;
        int hwm = -1;
        int lwm = -1;
        for (int d = difScoresSorted.length - 1; d >= 0; --d) {
            currentScore += difScoresSorted[d];
            if (hwm < 0 && currentScore >= target) {
                hwm = d;
            }
            
            if (lwm < 0 && hwm >= 0 && currentScore < target) {
                lwm = d;
            }
        }
        
        System.err.println("HWM=" + hwm);
        System.err.println("LWM=" + lwm);
        
        double expectedProp = 0.9;
        
        boolean[] flagArray = new boolean[negScores.length];
        for (int d = 0; d < flagArray.length; ++d) {
            flagArray[d] = posScores[d] > negScores[d];
        }
        double score = doScore(flagArray, negScores, posScores);
        System.err.println("Optimum score = " + score + " occupancy = " + occupancy(flagArray));
        int cycles = 0;
        while (cycles < 1000) {
          PICK_LOOP:
            do {
                if (Math.random() < 0.2) {
                    int dPlus = pick(flagArray, true);
                    int dMinus = pick(flagArray, false);
                    
                    double newScore = score - difScores[dPlus] + difScores[dMinus];
                    if (newScore > target) {
                        flagArray[dPlus] = false;
                        flagArray[dMinus] = true;
                        score = newScore;
                        break PICK_LOOP;
                    }
                } else {
	                int d = (int) (Math.random() * flagArray.length);
	                boolean current = flagArray[d];
	                if (current ^ (Math.random() < expectedProp)) {
	                    if (current) {
	                        double newScore = score - difScores[d];
	                        if (newScore > target) {
	                            flagArray[d] = false;
	                            score = newScore;
	                            break PICK_LOOP;
	                        }
	                    } else {
	                        double newScore = score + difScores[d];
	                        if (newScore > target) {
	                            flagArray[d] = true;
	                            score = newScore;
	                            break PICK_LOOP;
	                        }
	                    }
	                }
                }
            } while (true);
        	   ++cycles;
        }
        
        System.err.println("Got score " + score + " with occupancy " + occupancy(flagArray));
    }
    
    private static int pick(boolean[] flagArray, boolean val) {
        do {
            int d = (int) (Math.random() * flagArray.length);
            if (flagArray[d] == val) {
                return d;
            }
        } while (true);
    }
    
    private static double doScore(boolean[] flagArray, double[] negScore, double[] posScore) {
        double score = 0;
        for (int d = 0; d < flagArray.length; ++d) {
            score += (flagArray[d] ? posScore : negScore)[d];
        }
        return score;
    }
    
    private static double occupancy(boolean[] flagArray) {
        int num = 0;
        for (int d = 0; d < flagArray.length; ++d) {
            if (flagArray[d]) ++num;
        }
        return ((1.0 * num) / flagArray.length);
    }
}
