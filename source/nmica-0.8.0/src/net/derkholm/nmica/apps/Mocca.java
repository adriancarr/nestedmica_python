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

package net.derkholm.nmica.apps;

import net.derkholm.nmica.maths.DoubleFunction;
import net.derkholm.nmica.maths.IdentityDoubleFunction;
import net.derkholm.nmica.model.*;
import net.derkholm.nmica.model.logitseq.*;
import net.derkholm.nmica.model.motif.*;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;
import net.derkholm.nmica.trainer.EvaluationManager;
import net.derkholm.nmica.trainer.LocalEvaluationManager;
import net.derkholm.nmica.trainer.RandomDecopTrainer;
import net.derkholm.nmica.trainer.Trainer;
import net.derkholm.nmica.trainer.distributed.DistributedEvaluationManager;
import net.derkholm.nmica.utils.CliTools;
import net.derkholm.nmica.utils.ConfigurationException;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.io.FastaFormat;
import org.biojava.bio.seq.io.SeqIOTools;

import java.io.*;
import java.lang.reflect.Method;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Motif-oriented classification/clustering (MOCCA).
 * 
 * @author thomas
 */
public class Mocca {
    private int targetLength = 10;
    private int ensembleSize = -1;
    private int numMotifs = 10;
    private boolean revComp = false;
    private String snapshotFile = null;
    private String outFile = "classmotifs.xms";
    private int sampleInterval = 1000;
    private String checkpoint = null;
    private int keepCheckpoints = 2;
    private File restartFromCheckpoint = null;
    private int checkpointInterval = 10000;
    private int workerThreads = 1;
    private int maxCycles = 0;
    private SequenceDB[] seqDBs;
    private double crossOverProb = 0.3;
    private double simplexMaximumScale = 2.0;
    private double replaceComponentProb = 0.2;
    private double minClip = 0.002;
    private double maxClip = 1.0;
    private int minContributionMoves = -1;
    private int minContributionProposals = -1;
    private boolean distributed = false;
    private int port = 0;
    
    public void setCheckpoint(String checkpoint) {
        this.checkpoint = checkpoint;
    }
    public void setCheckpointInterval(int checkpointInterval) {
        this.checkpointInterval = checkpointInterval;
    }
    public void setCrossOverProb(double crossOverProb) {
        this.crossOverProb = crossOverProb;
    }
    public void setEnsembleSize(int ensembleSize) {
        this.ensembleSize = ensembleSize;
    }
    public void setKeepCheckpoints(int keepCheckpoints) {
        this.keepCheckpoints = keepCheckpoints;
    }
    public void setMaxCycles(int maxCycles) {
        this.maxCycles = maxCycles;
    }
    public void setMinContributionMoves(int minContributionMoves) {
        this.minContributionMoves = minContributionMoves;
    }
    public void setMinContributionProposals(int minContributionProposals) {
        this.minContributionProposals = minContributionProposals;
    }
    public void setNumMotifs(int numMotifs) {
        this.numMotifs = numMotifs;
    }
    public void setOutFile(String outFile) {
        this.outFile = outFile;
    }
    public void setReplaceComponentProb(double replaceComponentProb) {
        this.replaceComponentProb = replaceComponentProb;
    }
    public void setRestartFromCheckpoint(File restartFromCheckpoint) {
        this.restartFromCheckpoint = restartFromCheckpoint;
    }
    public void setRevComp(boolean revComp) {
        this.revComp = revComp;
    }
    public void setSampleInterval(int sampleInterval) {
        this.sampleInterval = sampleInterval;
    }
    public void setSeqs(File[] files)
    		throws Exception
    {
	    seqDBs = new SequenceDB[files.length];
	    for (int f = 0; f < files.length; ++f) {
	        seqDBs[f] = loadDB(files[f]);
	    }
	}
    private static SequenceDB loadDB(File f)
	    throws Exception
	{
        Pattern LABEL_PATTERN = Pattern.compile("label=(1|-1)");
        
	    SequenceIterator si = SeqIOTools.readFastaDNA(new BufferedReader(new FileReader(f)));
	    SequenceDB seqDB = new HashSequenceDB();
	    while (si.hasNext()) {
	        Sequence seq = si.nextSequence();
	        String descLine = (String) seq.getAnnotation().getProperty(FastaFormat.PROPERTY_DESCRIPTIONLINE);
	        Matcher matcher = LABEL_PATTERN.matcher(descLine);
	        if (matcher.find()) {
	            int label = Integer.parseInt(matcher.group(1));
	            seq.getAnnotation().setProperty("mocca.label", new Integer(label));
	        }
	        seqDB.addSequence(seq);
	    }
	    return seqDB;
	}
    public void setSeqDBs(SequenceDB[] seqDBs) {
        this.seqDBs = seqDBs;
    }
    public void setSimplexMaximumScale(double simplexMaximumScale) {
        this.simplexMaximumScale = simplexMaximumScale;
    }
    public void setSnapshotFile(String snapshotFile) {
        this.snapshotFile = snapshotFile;
    }
    public void setTargetLength(int targetLength) {
        this.targetLength = targetLength;
    }
    public void setWorkerThreads(int workerThreads) {
        this.workerThreads = workerThreads;
    }
    
    public static void main(String[] args) 
    		throws Exception
    {
        Mocca app = new Mocca();
        args = CliTools.configureBean(app, args);
        app.run(args);
    }
    
    public void run(String[] args)
    		throws Exception
    {
        Trainer trainer;
        
        boolean needsInit;
        if (restartFromCheckpoint != null) {
            ObjectInputStream ois = new ObjectInputStream(new FileInputStream(restartFromCheckpoint));
            trainer = (Trainer) ois.readObject();
            ois.close();
            
            numMotifs = trainer.getComponents();
            
            needsInit = false;
        } else {
            if (seqDBs == null) {
                throw new ConfigurationException("You must specify a -seqs option");
            }
            
            Datum[] data;
            {
                Set<String> allIds = new HashSet<String>();
                for (int s = 0; s < seqDBs.length; ++s) {
                    allIds.addAll((Set<String>) seqDBs[s].ids());
                }
                
                List<Datum> dl = new ArrayList<Datum>();
                for (String id: allIds) {
                    Object[] datumSeqs = new Object[seqDBs.length];
                    for (int s = 0; s < seqDBs.length; ++s) {
                        if (seqDBs[s].ids().contains(id)) {
                            datumSeqs[s] = seqDBs[s].getSequence(id);
                        }
                    }
                    dl.add(new SimpleDatum(id, datumSeqs));
                }
                data = dl.toArray(new Datum[dl.size()]);
            }
            
            
            // sane defaults for some stuff
            
            if (ensembleSize < 0) {
                ensembleSize = Math.max(200, (int) (1000.0 / numMotifs));
            }
            if (minContributionMoves < 0) {
                minContributionMoves = numMotifs * 2;
            }
            if (minContributionProposals < 0) {
                minContributionProposals = numMotifs * 8;
            }
            
            ContributionPrior prior;
            WeightMatrixPrior wmPrior;
            wmPrior = new MotifClippedSimplexPrior(org.biojava.bio.seq.DNATools.getDNA(), targetLength, targetLength, targetLength, minClip, maxClip);
            prior = new WeightedWeightMatrixPrior(wmPrior);
                
            ContributionSampler sampler;
            {
                MultiplexContributionSampler mcsSampler = new MultiplexContributionSampler();
                mcsSampler.addSampler(new WeightedWeightMatrixMatrixSampler(
                        new SymbolMassScalingSampler(simplexMaximumScale)
                ), 16.0);
                /*
                mcsSampler.addSampler(new WeightedWeightMatrixMatrixSampler(
                        new SpinSampler(wmPrior)
                ), 2.0);
                */
                /*
                mcsSampler.addSampler(new WeightedWeightMatrixMatrixSampler(
                        new ZapSampler(wmPrior)
                ), 4.0);
                mcsSampler.addSampler(new WeightedWeightMatrixMatrixSampler(
                        new IndelSampler(wmPrior)
                ), 2.0);
                */
                mcsSampler.addSampler(new WeightedWeightMatrixWeightSampler(), 5.0);
                sampler = mcsSampler;
            }
            
            MixPolicy mixPolicy = new FlatMixPolicy();
            DoubleFunction mixTransferFunction = IdentityDoubleFunction.INSTANCE;
            
            Facette[] facettes = new Facette[seqDBs.length];
            for (int f = 0; f < seqDBs.length; ++f) {
	            facettes[f] = new LogisticSequenceFacette();
            }
            ContributionGroup cg = new SimpleContributionGroup("motifs", WeightedWeightMatrix.class);
            SimpleFacetteMap facetteMap = new SimpleFacetteMap(
                new ContributionGroup[] {cg},
                facettes
            );
            for (int f = 0; f < facettes.length; ++f) {
                facetteMap.setContributesToFacette(cg, facettes[f], true);
            }

	        trainer = new RandomDecopTrainer(
	                facetteMap,
	                data,
	                numMotifs,
	                new ContributionPrior[] { prior },
	                new ContributionSampler[] { sampler },
	                mixPolicy,
	                ensembleSize
	        );
	        ((RandomDecopTrainer) trainer).setMinMixtureMoves(0);
	        ((RandomDecopTrainer) trainer).setMinContributionMoves(minContributionMoves);
	        ((RandomDecopTrainer) trainer).setMinMixtureProposals(0);
	        ((RandomDecopTrainer) trainer).setMinContributionProposals(minContributionProposals);
	        trainer.setIgnoreMixturePrior(true);
	            
            trainer.setCrossOverProb(crossOverProb);
            trainer.setReplaceComponentProb(replaceComponentProb);
            needsInit = true;
        }
        
        EvaluationManager em;
        if (distributed) {
            em = new DistributedEvaluationManager(port, false, false);
        } else {
            em = new LocalEvaluationManager(workerThreads);
        }
        trainer.setEvaluationManager(em);
        
        if (needsInit) {
            trainer.init();
        }
        
        {
            try {
                Method shortCircuitCall = org.biojava.utils.ChangeSupport.class.getMethod(
                    "setGlobalChangeBypass", 
                    new Class[] {Boolean.TYPE}
                );
                shortCircuitCall.invoke(null, new Object[] {Boolean.TRUE});
            } catch (NoSuchMethodException ex) {
                System.err.println("Short-circuiting isn't available");
            }
        }
        
        List<File> checkpointFiles = new ArrayList<File>();
        
        long lastCycleTime = 0;
        while (true) {
            Trainer.WeightedModel twm = trainer.next();
            int c = trainer.getCycle();
            if (snapshotFile != null && (c % sampleInterval == 0)) {
                writeMotifs(trainer.getBestModel(), snapshotFile + '.' + c + ".xms");
            }
            long currentCycleTime = System.currentTimeMillis();
            System.out.print("" + c + '\t' + twm.getWeight() + '\t' + twm.getLikelihood());
            System.out.print("\t" + trainer.getLikelihoodIQR());
            System.out.print("\t" + (currentCycleTime - lastCycleTime));
            lastCycleTime = currentCycleTime;
            System.out.println();
            
            if (checkpoint != null && (c % checkpointInterval) == 0) {
                File cpFile = new File(checkpoint + '.' + c + ".jos");
                ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(cpFile));
                oos.writeObject(trainer);
                oos.close();
                checkpointFiles.add(cpFile);
                while (checkpointFiles.size() > keepCheckpoints) {
                    File oldcp = checkpointFiles.remove(0);
                    oldcp.delete();
                }
            }
            
            if ((maxCycles > 0 && c >= maxCycles) || trainer.getLikelihoodIQR() < 0.01) {
                if (outFile != null) {
                    writeMotifs(trainer.getBestModel(), outFile);
                }
                
                if (em instanceof DistributedEvaluationManager) {
                    ((DistributedEvaluationManager) em).shutdown();
                }
                return;
            }
        } 
    }
    
    private void writeMotifs(MultiICAModel model, String fileName)
		throws Exception
	{
		ContributionGroup motifCG = model.getFacetteMap().getContributionGroups()[0];
		List<Motif> motifList = new ArrayList<Motif>();
		for (int m = 0; m < model.getComponents(); ++m) {
		    WeightedWeightMatrix wwm = (WeightedWeightMatrix) model.getContribution(motifCG, m).getItem();
		    
		    Motif motif = new Motif();
		    motif.setName("motif" + m);
		    motif.setWeightMatrix(wwm.getWeightMatrix());
		    motif.getAnnotation().setProperty("mocca.weight", new Double(wwm.getWeight()));
		    motifList.add(motif);
		}
		Motif[] motifs = motifList.toArray(new Motif[0]);
		
		FileOutputStream fos = new FileOutputStream(fileName);
		MotifIOTools.writeMotifSetXML(fos, motifs);
		fos.close();
	}
}
