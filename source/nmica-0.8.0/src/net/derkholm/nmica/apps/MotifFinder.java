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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.lang.reflect.Method;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TreeSet;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamReader;
import javax.xml.stream.XMLStreamWriter;

import net.derkholm.nmica.build.NMApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.maths.DoubleFunction;
import net.derkholm.nmica.maths.IdentityDoubleFunction;
import net.derkholm.nmica.maths.LogisticFunction;
import net.derkholm.nmica.maths.NativeMath;
import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.matrix.SimpleMatrix2D;
import net.derkholm.nmica.model.BinaryMixPolicy;
import net.derkholm.nmica.model.ConstantMixPolicy;
import net.derkholm.nmica.model.ContributionGroup;
import net.derkholm.nmica.model.ContributionItem;
import net.derkholm.nmica.model.ContributionPrior;
import net.derkholm.nmica.model.ContributionSampler;
import net.derkholm.nmica.model.Datum;
import net.derkholm.nmica.model.FlatMixPolicy;
import net.derkholm.nmica.model.HistoryThread;
import net.derkholm.nmica.model.MixPolicy;
import net.derkholm.nmica.model.MultiICAModel;
import net.derkholm.nmica.model.MultiplexContributionSampler;
import net.derkholm.nmica.model.OneOfManyMixPolicy;
import net.derkholm.nmica.model.RealMixPolicy;
import net.derkholm.nmica.model.SimpleContributionGroup;
import net.derkholm.nmica.model.SimpleContributionItem;
import net.derkholm.nmica.model.SimpleDatum;
import net.derkholm.nmica.model.SimpleFacetteMap;
import net.derkholm.nmica.model.WeightedColumnMixPolicy;
import net.derkholm.nmica.model.motif.Mosaic;
import net.derkholm.nmica.model.motif.MosaicIO;
import net.derkholm.nmica.model.motif.MosaicSequenceBackground;
import net.derkholm.nmica.model.motif.MosaicTools;
import net.derkholm.nmica.model.motif.MotifClippedSimplexPrior;
import net.derkholm.nmica.model.motif.MotifFacette;
import net.derkholm.nmica.model.motif.RetrimSampler;
import net.derkholm.nmica.model.motif.SequenceBackground;
import net.derkholm.nmica.model.motif.SlideSampler;
import net.derkholm.nmica.model.motif.SymbolMassScalingSampler;
import net.derkholm.nmica.model.motif.WeightMatrixPrior;
import net.derkholm.nmica.model.motif.ZapSampler;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;
import net.derkholm.nmica.trainer.EvaluationManager;
import net.derkholm.nmica.trainer.LocalEvaluationManager;
import net.derkholm.nmica.trainer.MixtureResamplingTrainer;
import net.derkholm.nmica.trainer.QueuedDecopTrainer;
import net.derkholm.nmica.trainer.RandomDecopTrainer;
import net.derkholm.nmica.trainer.Trainer;
import net.derkholm.nmica.trainer.distributed.DistributedEvaluationManager;
import net.derkholm.nmica.utils.ArrayTools;
import net.derkholm.nmica.utils.CollectTools;

import org.biojava.bio.Annotation;
import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.io.FastaFormat;
import org.biojava.bio.seq.io.SeqIOTools;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.ConfigurationException;
import org.bjv2.util.cli.Option;
import org.bjv2.util.cli.UserLevel;

@App(overview="The NestedMICA motif finder.  NestedMICA uses nested sampling to learn a dictionary of over-represented sequence motifs in a set of DNA or Protein sequence data.", generateStub=true)
@NMApp(launchName="nminfer", vm=VirtualMachine.SERVER)
public class MotifFinder {
    private boolean counted = false;
	
    private FiniteAlphabet alphabet = DNATools.getDNA(); 
    private String mixtureType = "binary";
    private int maxLength = 10;
    private int minLength = -1;
    private int extraLength = 2;
    private int ensembleSize = -1;
    private int numMotifs = 1;
    private int minSeqLength = 0;
    private boolean revComp = false;
    private SequenceBackground[] backgroundModels;
    private String sampleFile = null;
    private String outFile = "motifs.xms";
    private int sampleInterval = 1000;
    private String checkpoint = null;
    private int keepCheckpoints = 2;
    private File restartFromCheckpoint = null;
    private int checkpointInterval = 10000;
    private int workerThreads = 1;
    private int maxCycles = 0;
    private File[] seqDBfiles;
    private SequenceDB[] seqDBs;
    private double crossOverProb = 0.3;
    private double uncountedExpectation = 1.0;
    private double simplexMaximumScale = 2.0;
    private double replaceComponentProb = 0.2;
    private double edgePrune = 0.0;
    private boolean discriminate = false;
    private boolean restrictCrossover = true;
    
    private File writeBackground = null;
    private int backgroundClasses = 1;
    private int backgroundOrder = 0;
	    
    private WeightMatrix[] knownMotifs = new WeightMatrix[0];
    
    private String mixtureUpdate = "resample";
    
    private int minMixtureMoves = 2;
    private int minContributionMoves = -1;
    private int minMixtureProposals = -1;
    private int minContributionProposals = -1;
    private int mixtureDecopSessions = -1;
    private double mixtureDecopFraction = 0.8;
    private boolean ignoreMixturePrior = false;
    
    private boolean cluster = false;
    private double clusterIn = 0.01;
    private double clusterOut = 0.02;
    
    private double minClip = Double.NEGATIVE_INFINITY;
    private double maxClip = 1.0;
    
    private double expectedUsageFraction = 0.5;
    private double mixStepSize = 0.2;
    
    private boolean distributed = false;
    private int port = 0;
    private boolean crabDebug = false;
    private boolean throughtputMonitor = false;
    private double crabSignifier = -1;
    private double crabRate = -1;
    private boolean disableBalancer = true;
    
    private boolean cardinalityOut = false;
    private int logInterval = 10;

    private boolean terminateOnConvergance = true;

    @Option(help="Write out a background model file (only applicable if -backgroundModel not specifified)", optional=true)
    public void setWriteBackground(File f) {
    	this.writeBackground = f;
    }
    
    @Option(help="Number of classes to use for background models (only used if -backgroundModel not specified, default=1)", optional=true)
    public void setBackgroundClasses(int i) {
    	this.backgroundClasses = i;
    }
    
    @Option(help="Order of background models (only used if -backgroundModel not specified, default=0)", optional=true)
    public void setBackgroundOrder(int i) {
    	this.backgroundOrder = i;
    }
    
    @Option(help="Disable the load-balancer", userLevel=UserLevel.DEBUG, optional=true, negation="enableBalancer")
    public void setDisableBalancer(boolean b) {
    	this.disableBalancer = b;
    }
    
    @Option(help="Step size when adjusting CRAB weights", userLevel=UserLevel.DEBUG, optional=true)
    public void setCrabRate(double i) {
    	this.crabRate = i;
    }
    
    @Option(help="Smallest workload to consider when adjusting CRAB weights", userLevel=UserLevel.DEBUG, optional=true)
    public void setCrabSignifier(double i) {
    	this.crabSignifier = i;
    }
    
    @Option(help="Enable monitoring of load-balancer weights (default=false)", userLevel=UserLevel.DEBUG, optional=true)
    public void setCrabMonitor(boolean b) {
    	this.crabDebug = true;
    }
    
    @Option(help="Enable throughput instrumentation in distributed mode (default=false)", userLevel=UserLevel.DEBUG, optional=true)
    public void setThroughputMonitor(boolean b) {
    	this.throughtputMonitor = b;
    }
    
    @Option(help="Prevent crossovers which might create duplicate motifs (default=true)", userLevel=UserLevel.EXPERT, optional=true)
    public void setRestrictCrossover(boolean b) {
    	this.restrictCrossover = b;
    }
    
    @Option(help="Counted motif DP (slow!)", userLevel=UserLevel.EXPERT, optional=true)
    public void setCounted(boolean b) {
    	this.counted = b;
    }
    
    @Option(help="Discriminitive likelihood calculations (test purposes only)", userLevel=UserLevel.EXPERT, optional=true)
    public void setDiscriminate(boolean b) {
    	this.discriminate = b;
    }
    
    @Option(help="A set of known motifs.  These are softmasked when learning new motifs.", userLevel=UserLevel.EXPERT, optional=true)
    public void setKnownMotifs(InputStream is)
        throws Exception
    {
        Motif[] km = MotifIOTools.loadMotifSetXML(is);
        knownMotifs = new WeightMatrix[km.length];
        for (int i = 0; i < km.length; ++i) {
            knownMotifs[i] = km[i].getWeightMatrix();
        }
    }

    @Option(help="Alphabet of the input sequences (default=DNA)", optional=true)
    public void setAlphabet(MicaAlphabet alpha) 
    {
    	this.alphabet = alpha.biojavaAlphabet();
    }
    
    @Option(help="Automatic termination (default=true)", optional=true, userLevel=UserLevel.EXPERT)
    public void setTerminateOnConvergance(boolean b) {
        this.terminateOnConvergance = b;
    }
    public void setCardinalityOut(boolean b) {
        this.cardinalityOut = b;
    }
    @Option(help="The number of cycles between likelihood log points.", optional=true)
    public void setLogInterval(int i) {
        this.logInterval = i;
    }
    @Option(help="The strategy for updating the ICA mixing matrix.", optional=true, userLevel=UserLevel.EXPERT)
    public void setMixtureUpdate(String b) {
        this.mixtureUpdate = b;
    }
    
    @Option(help="Don't calculate the prior over the mixing matrix", optional=true, userLevel=UserLevel.EXPERT)
    public void setIgnoreMixturePrior(boolean b) {
        this.ignoreMixturePrior = b;
    }
    
    @Option(help="The number of steps to take when resampling the ICA mixing matrix.", optional=true, userLevel=UserLevel.EXPERT)
    public void setMixtureDecopSessions(int i) {
        this.mixtureDecopSessions = i;
    }

    public void setMixtureDecopFraction(double d) {
        this.mixtureDecopFraction = d;
    }

    public void setKeepCheckpoints(int i) {
        this.keepCheckpoints = i;
    }
    
    @Option(help="The expected fraction of motifs which are used in each sequence", optional=true, userLevel=UserLevel.EXPERT)
    public void setExpectedUsageFraction(double d) {
        this.expectedUsageFraction = d;
    }
    
    @Option(help="Operate in distributed mode.  Run dlepnode processes to perform calculations on behalf of this motiffinder process", optional=true)
    public void setDistributed(boolean b) {
        this.distributed = b;
    }
    
    @Option(help="The network port to use in distributed mode", optional=true)
    public void setPort(int p) {
        this.port = p;
    }
    
    public void setMixStepSize(double d) {
        this.mixStepSize = d;
    }
    
    @Option(help="The type of ICA mixing matrix to use (binary|flat)", optional=true, userLevel=UserLevel.EXPERT)
    public void setMixtureType(String s) {
        this.mixtureType = s;
    }
    
    @Option(help="Minimum weight for a symbol", userLevel = UserLevel.EXPERT, optional=true)
    public void setMinClip(double d) {
        this.minClip = d;
    }
    
    public void setMaxClip(double d) {
        this.maxClip = d;
    }

    @Option(help="Prefer motifs which occur in clusters on the sequences", userLevel=UserLevel.EXPERT, optional=true)
	public void setCluster(boolean cluster) {
		this.cluster = cluster;
	}
	
	public void setClusterIn(double clusterIn) {
		this.clusterIn = clusterIn;
	}
	public void setClusterOut(double clusterOut) {
		this.clusterOut = clusterOut;
	}
    
    public void setEdgePrune(double d) {
        this.edgePrune = d;
    }
    
    public void setReplaceComponentProb(double d) {
        this.replaceComponentProb = d;
    }
    
    public void setSimplexMaximumScale(double d) {
        this.simplexMaximumScale = d;
    }
    
    @Option(help="The expected number of motif occurrences per sequence", userLevel=UserLevel.EXPERT, optional=true)
    public void setUncountedExpectation(double i) {
    		this.uncountedExpectation = i;
    }
    
    @Option(help="...", userLevel=UserLevel.DEBUG, optional=true)
    public void setCrossOverProb(double d) {
        this.crossOverProb = d;
    }

    @Option(help="The sequences to analyse", optional=true)
    public void setSeqs(File[] files)
        throws Exception
    {
    	this.seqDBfiles = files;
    }
    
    @Option(help="Maximum number of cycles to run", optional=true)
    public void setMaxCycles(int i) {
        this.maxCycles = i;
    }
    
    @Option(help="Number of threads to use.  For best performance, set this equal to the number of CPUs in your computer", optional=true)
    public void setThreads(int i) {
        this.workerThreads = i;
    }
    
    @Option(help="Interval (in cycles) between sampling the state of the motif finder", optional=true)
    public void setSampleInterval(int i) {
        this.sampleInterval = i;
    }
    
    @Option(help="Restart a previously saved checkpoint file.  Most other options are silently ignored when this is used", optional=true)
    public void setRestartFromCheckpoint(File f) {
        this.restartFromCheckpoint = f;
    }
    
    @Option(help="Base name for checkpoint files", optional=true)
    public void setCheckpoint(String f) {
        this.checkpoint = f;
    }
    
    @Option(help="Interval (in cycles) between writing checkpoint files", optional=true)
    public void setCheckpointInterval(int interval) {
        this.checkpointInterval = interval;
    }
    
    @Option(help="Allow motifs to occur in either orientation", optional=true)
    public void setRevComp(boolean b) {
        this.revComp = b;
    }
    
    @Option(help="Base name of files used to periodically sample the motif-finder's current state", optional=true)
    public void setSampleFile(String s) {
        this.sampleFile = s;
    }
    
    @Option(help="Name of file for recording the final set of motifs", optional=true)
    public void setOut(String s) {
        this.outFile = s;
    }
    
    @Option(help="The background model to use.  You can create new background model files with the makemosaicbg program", optional=true)
    public void setBackgroundModel(InputStream[] is)
        throws Exception
    {
        backgroundModels = new SequenceBackground[is.length];
        XMLInputFactory factory = XMLInputFactory.newInstance();
        
        for (int b = 0; b < is.length; ++b) {
        	try {
        		XMLStreamReader r = factory.createXMLStreamReader(is[b]);
        		Mosaic m = MosaicIO.readMosaic(r);
        		backgroundModels[b] = new MosaicSequenceBackground(m.getDistributions(), m.getTransition());
        		r.close();
        	} catch (Exception ex) {
        		// ex.printStackTrace();
        		throw new Exception("Error loading background model.  If you are using a background model created with an earlier version of NestedMICA, please try using the nmconvertbg program", ex);
        	}
        }
    }
    
    @Option(help="...", userLevel=UserLevel.EXPERT, optional=true)
    public void setMinMixtureMoves(int m) {
        this.minMixtureMoves = m;
    }
    
    @Option(help="...", userLevel=UserLevel.EXPERT, optional=true)
    public void setMinContributionMoves(int m) {
        this.minContributionMoves = m;
    }
    
    /**
     * @param minContributionProposals The minContributionProposals to set.
     */
    @Option(help="...", userLevel=UserLevel.EXPERT, optional=true)
    public void setMinContributionProposals(int minContributionProposals) {
        this.minContributionProposals = minContributionProposals;
    }
    /**
     * @param minMixtureProposals The minMixtureProposals to set.
     */
    @Option(help="...", userLevel=UserLevel.EXPERT, optional=true)
    public void setMinMixtureProposals(int minMixtureProposals) {
        this.minMixtureProposals = minMixtureProposals;
    }
    
    @Option(help="The maximum length of learned motifs", optional=true)
    public void setMaxLength(int i) {
        this.maxLength = i;
    }
    
    @Option(help="The minimum length of learned motifs (if not specified, this is equal to maxLength)", optional=true)
    public void setMinLength(int i) {
        this.minLength = i;
    }
    
    @Option(help="Extra flanking bases which support the calculation in variable-length mode", userLevel=UserLevel.EXPERT, optional=true)
    public void setExtraLength(int i) {
    	this.extraLength = i;
    }
    
    @Option(help="Size of the nested sampling ensemble", userLevel=UserLevel.EXPERT, optional=true)
    public void setEnsembleSize(int i) {
        this.ensembleSize = i;
    }
    
    @Option(help="The number of motifs to learn", optional=true)
    public void setNumMotifs(int i) {
        this.numMotifs = i;
    }

    @Option(help="The minimum length for sequences used during training", optional=true)
    public void setMinSeqLength(int i) {
        this.minSeqLength = i;
    }
    
    private double scale(int n, int x0, double y0, int x1, double y1) {
    	// n = n*n;
    	// x0 = x0 * x0;
    	// x1 = x1 * x1;
    	if (n <= x0) {
    		return y0;
    	} else if (n >= x1) {
    		return y1;
    	} else {
    		return y0 + ((1.0 * n) - x0) / (x1 - x0) * (y1 - y0);
    	}
    }
    
    public void main(String[] args)
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
        	
            if (seqDBfiles == null) {
                throw new ConfigurationException("You must specify a -seqs option");
            }
        	
        	seqDBs = new SequenceDB[seqDBfiles.length];
            for (int f = 0; f < seqDBfiles.length; ++f) {
                seqDBs[f] = loadDB(seqDBfiles[f], alphabet);
            }
  
            if (backgroundModels == null) {
                System.err.println("*** Warning: no background model specified, generating a background model using the input sequences");
                backgroundModels = new SequenceBackground[seqDBs.length];
                double mosaicTransition = 0.005;
                for (int s = 0; s < seqDBs.length; ++s) {
                    Distribution[] patches = MosaicTools.optimizePatches(
                        seqDBs[s],
                        backgroundClasses,
                        backgroundOrder + 1, 
                        mosaicTransition,
                        alphabet == ProteinTools.getAlphabet()
                    );
                	backgroundModels[s] = new MosaicSequenceBackground(patches, mosaicTransition);
                	
                	if (writeBackground != null && s == 0) {
                        Mosaic mosaic = new Mosaic(patches, mosaicTransition);
                        mosaic.getAnnotation().setProperty("creator.name", "nminfer");
                        mosaic.getAnnotation().setProperty("creator.version", Version.VERSION);
                        mosaic.getAnnotation().setProperty("input", seqDBfiles[0].getName());
                        mosaic.getAnnotation().setProperty("date", new SimpleDateFormat().format(new Date()));
                        
                        XMLOutputFactory factory = XMLOutputFactory.newInstance();
                        XMLStreamWriter xw = factory.createXMLStreamWriter(new FileWriter(writeBackground));
                        MosaicIO.writeMosaic(xw, mosaic);
                        xw.close();
                	}
                }
            } else if (backgroundModels.length != seqDBs.length && backgroundModels.length != 1) {
                throw new ConfigurationException("Mismatch between number of -seqs and -backgroundModel arguments");
            }
            
            Datum[] data;
            {
                Set<String> allIds = new HashSet<String>();
                for (int s = 0; s < seqDBs.length; ++s) {
                    allIds.addAll((Set<? extends String>) seqDBs[s].ids());
                }
                
                List<Datum> dl = new ArrayList<Datum>();
                for (Iterator<String> i = allIds.iterator(); i.hasNext(); ) {
                    String id = i.next();
                    Object[] datumSeqs = new Object[seqDBs.length];
                    boolean gotDatumSeq = false;
                    for (int s = 0; s < seqDBs.length; ++s) {
                        if (seqDBs[s].ids().contains(id)) {
                            Sequence seq = seqDBs[s].getSequence(id);
                            
                            int ssOrder = 0;
                            SequenceBackground ssBg = backgroundModels.length > 1 ? backgroundModels[s] : backgroundModels[0];
                            if (ssBg instanceof MosaicSequenceBackground) {
                            	ssOrder = ((MosaicSequenceBackground) ssBg).getBackgroundDistributions()[0].getAlphabet().getAlphabets().size() - 1;
                            }
                            
                            if (seq.length() < maxLength + ssOrder || seq.length() < minSeqLength) {
                                System.err.println("*** Warning: discarding tiny sequence " + id);
                            } else {
                                datumSeqs[s] = seqDBs[s].getSequence(id);
                                gotDatumSeq = true;
                            }
                        }
                    }
                    if (gotDatumSeq) {
                        dl.add(new SimpleDatum(id, datumSeqs));
                    }
                }
                data = dl.toArray(new Datum[dl.size()]);
            }
            
            numMotifs += knownMotifs.length;
            
            // sane defaults for some stuff
            
            if (minLength < 0) {
            	minLength = maxLength;
            	System.err.printf("*** Warning: minLength not specified, assuming minLength=maxLength=%d%n", maxLength);
            } else if (minLength > maxLength) {
            	System.err.println("minLength must be less than or equal to maxLength");
            	return;
            }
            
            if (minClip < 0) {
            	if (alphabet == ProteinTools.getAlphabet()) {
            		minClip = 0.0000001;
            	} else {
            		minClip = 0.002;
            	}
            }
            if (ensembleSize < 0) {
            	if (alphabet == ProteinTools.getAlphabet()) {
            		ensembleSize = Math.max(500, (int) (1500.0 / numMotifs));
            	} else {
            		ensembleSize = Math.max(200, (int) (1000.0 / numMotifs));
            	}
            }
            if (minContributionMoves < 0) {
                minContributionMoves = (int) Math.ceil(numMotifs * scale(numMotifs, 1, 2, 10, 1));
            }
            if (minContributionProposals < 0) {
                minContributionProposals = (int) Math.ceil(numMotifs * scale(numMotifs, 1, 8, 10, 3));
            }
            if (minMixtureMoves < 0) {
                minMixtureMoves = (int) Math.ceil(0.01 * data.length * numMotifs);
            }
            if (minMixtureProposals < 0) {
                minMixtureProposals = (int) Math.ceil(0.25 * data.length * numMotifs);
            }
            if (mixtureDecopSessions < 0) {
                mixtureDecopSessions = (int) Math.ceil(0.4 * numMotifs);
            }
            if ("flat".equals(mixtureType)) {
                mixtureUpdate = "random";
                minMixtureProposals = 0;
                minMixtureMoves = 0;
            }
            
            ContributionPrior prior;
            WeightMatrixPrior wmPrior;
            prior = wmPrior = new MotifClippedSimplexPrior(alphabet, minLength, maxLength, maxLength + extraLength, minClip, maxClip);
                
            ContributionSampler sampler;
            {
                MultiplexContributionSampler mcsSampler = new MultiplexContributionSampler();
                mcsSampler.addSampler(new SymbolMassScalingSampler(simplexMaximumScale), 16.0);
                mcsSampler.addSampler(new SlideSampler(), 2.0);
                mcsSampler.addSampler(new ZapSampler(wmPrior), 4.0);
                // mcsSampler.addSampler(new IndelSampler(wmPrior), 2.0);  would be nice to replace this...
                if (minLength < maxLength) {
                	mcsSampler.addSampler(new RetrimSampler(wmPrior));
                }
                sampler = mcsSampler;
            }
            
            MixPolicy mixPolicy;
            DoubleFunction mixTransferFunction = IdentityDoubleFunction.INSTANCE;
            
            if (discriminate) {
            	Set<String> motifTagSet = new TreeSet<String>();
            	for (Datum d : data) {
            		Sequence s = (Sequence) d.getFacettedData()[0];
            		if (s != null) {
	            		String fdl = s.getAnnotation().getProperty(FastaFormat.PROPERTY_DESCRIPTIONLINE).toString();
	            		String[] toks = fdl.split(" ");
	            		for (int t = 1; t < toks.length; ++t) {
	            			motifTagSet.add(toks[t]);
	            		}
            		}
            	}
            	
            	String[] motifTags = motifTagSet.toArray(new String[0]);
            	numMotifs = motifTags.length;
            	if (numMotifs == 0) {
            		System.err.println("*** Couldn't find any motif tags in discriminitive dataset");
            		return;
            	} else {
            		System.err.println("*** In discriminitive mode, setting numMotifs=" + numMotifs);
            	}
            	
            	Matrix2D mm = new SimpleMatrix2D(data.length, numMotifs);
            	for (int di = 0; di < data.length; ++di) {
            		Sequence s = (Sequence) data[di].getFacettedData()[0];
            		if (s != null) {
	            		String fdl = s.getAnnotation().getProperty(FastaFormat.PROPERTY_DESCRIPTIONLINE).toString();
	            		String[] toks = fdl.split(" ");
	            		for (int t = 1; t < toks.length; ++t) {
	            			int mti = ArrayTools.indexOf(motifTags, toks[t]);
	            			mm.set(di, mti, 1.0);
	            		}
            		}
            	}
            	mixPolicy = new ConstantMixPolicy(mm);
            	
            	mixtureUpdate = "resample";
            	mixtureDecopSessions = 0;
            	crossOverProb = 0;
            } else {
	            if ("oneOfMany".equals(mixtureType)) {
	                mixPolicy = new OneOfManyMixPolicy();
	            } else if ("flat".equals(mixtureType)) {
	                mixPolicy = new FlatMixPolicy();
	            } else if ("binary".equals(mixtureType)) {
	                mixPolicy = new BinaryMixPolicy(expectedUsageFraction);
	            } else if ("logit".equals(mixtureType)) {
	                mixPolicy = new RealMixPolicy(-20, 20, mixStepSize);
	                mixTransferFunction = LogisticFunction.INSTANCE;
	            } else if (mixtureType.startsWith("weighted:")) {
	                List<Double> weights = new ArrayList<Double>();
	                for (StringTokenizer toke = new StringTokenizer(mixtureType.substring(9), ",="); toke.hasMoreTokens(); ) {
	                    int count = Integer.parseInt(toke.nextToken());
	                    Double weight = new Double(toke.nextToken());
	                    for (int c = 0; c < count; ++c) {
	                        weights.add(weight);
	                    }
	                }
	                mixPolicy = new WeightedColumnMixPolicy(CollectTools.toDoubleArray(weights));
	                if (numMotifs != weights.size()) {
	                    System.err.println("*** Warning: in weighted binary mixture mode but weight vector doesn't match numMotifs");
	                    System.err.println("*** Setting numMotifs=" + weights.size());
	                    numMotifs = weights.size();
	                }
	            } else {
	                System.err.println("Unsupported mixture type: " + mixtureType);
	                return;
	            }
            }
            
            MotifFacette[] facettes = new MotifFacette[seqDBs.length];
            for (int f = 0; f < seqDBs.length; ++f) {
	                facettes[f] = new MotifFacette(
	                    backgroundModels[backgroundModels.length > 1 ? f : 0],
	                    0.0,
	                    true,
	                    !counted,
	                    revComp,
	                    uncountedExpectation,
	                    false,
	                    edgePrune,
	                    alphabet
	                );
	                facettes[f].setCluster(cluster);
	                facettes[f].setClusterIn(clusterIn);
	                facettes[f].setClusterOut(clusterOut);
	                facettes[f].setMixTransferFunction(mixTransferFunction);
	                facettes[f].setDiscriminate(discriminate);
 
            }
            ContributionGroup cg = new SimpleContributionGroup("motifs", WeightMatrix.class);
            SimpleFacetteMap facetteMap = new SimpleFacetteMap(
                new ContributionGroup[] {cg},
                facettes
            );
            for (int f = 0; f < facettes.length; ++f) {
                facetteMap.setContributesToFacette(cg, facettes[f], true);
            }
            
            if ("max".equals(mixtureUpdate)) {
                MixtureResamplingTrainer mrt = new MixtureResamplingTrainer(
                        facetteMap,
        	                data,
        	                numMotifs,
        	                new ContributionPrior[] { prior },
        	                new ContributionSampler[] { sampler },
        	                mixPolicy,
        	                ensembleSize
        	           );
                   mrt.setMinContributionMoves(minContributionMoves);
                   mrt.setMinContributionProposals(minContributionProposals);
                   mrt.setMixtureDecopSessions(mixtureDecopSessions);
                   mrt.setMoveFraction(1.0);
                   mrt.setProposalFraction(0.01);
                   trainer = mrt;
            } else if ("weakResample".equals(mixtureUpdate)) {
                MixtureResamplingTrainer mrt = new MixtureResamplingTrainer(
                        facetteMap,
        	                data,
        	                numMotifs,
        	                new ContributionPrior[] { prior },
        	                new ContributionSampler[] { sampler },
        	                mixPolicy,
        	                ensembleSize
        	           );
                   mrt.setMinContributionMoves(minContributionMoves);
                   mrt.setMinContributionProposals(minContributionProposals);
                   mrt.setMixtureDecopSessions(mixtureDecopSessions);
                   mrt.setMoveFraction(0.02);
                   mrt.setProposalFraction(0.1);
                   mrt.setIgnoreMixturePrior(ignoreMixturePrior);
                   trainer = mrt;
            } else if ("resample".equals(mixtureUpdate)) {
                MixtureResamplingTrainer mrt = new MixtureResamplingTrainer(
                    facetteMap,
    	                data,
    	                numMotifs,
    	                new ContributionPrior[] { prior },
    	                new ContributionSampler[] { sampler },
    	                mixPolicy,
    	                ensembleSize
    	           );
               mrt.setMinContributionMoves(minContributionMoves);
               mrt.setMinContributionProposals(minContributionProposals);
               mrt.setMixtureDecopSessions(mixtureDecopSessions);
               mrt.setIgnoreMixturePrior(ignoreMixturePrior);
               
               if (discriminate) {
            	   mrt.setPermuteProbability(0);
               }
               
               trainer = mrt;
            } else if ("queue".equals(mixtureUpdate)) {
                QueuedDecopTrainer qdt = new QueuedDecopTrainer(
                    facetteMap,
    	                data,
    	                numMotifs,
    	                new ContributionPrior[] { prior },
    	                new ContributionSampler[] { sampler },
    	                mixPolicy,
    	                ensembleSize
    	           );
                qdt.setMinContributionMoves(minContributionMoves);
                qdt.setMinContributionProposals(minContributionProposals);
                qdt.setMixtureDecopSessions(mixtureDecopSessions);
                qdt.setMixtureFractionPerSession(mixtureDecopFraction);
                qdt.setIgnoreMixturePrior(ignoreMixturePrior);
                trainer = qdt;
            } else if ("random".equals(mixtureUpdate)){
	            trainer = new RandomDecopTrainer(
	                facetteMap,
	                data,
	                numMotifs,
	                new ContributionPrior[] { prior },
	                new ContributionSampler[] { sampler },
	                mixPolicy,
	                ensembleSize
	            );
	            ((RandomDecopTrainer) trainer).setMinMixtureMoves(minMixtureMoves);
	            ((RandomDecopTrainer) trainer).setMinContributionMoves(minContributionMoves);
	            ((RandomDecopTrainer) trainer).setMinMixtureProposals(minMixtureProposals);
	            ((RandomDecopTrainer) trainer).setMinContributionProposals(minContributionProposals);
	            trainer.setIgnoreMixturePrior(ignoreMixturePrior);
            } else {
                System.err.println("Unknown mixtureUpdate type: " + mixtureUpdate);
                return;
            }
            
            if (knownMotifs.length > 0) {
                ContributionItem[] seedCis = new ContributionItem[knownMotifs.length];
                for (int i = 0; i < knownMotifs.length; ++i) {
                    seedCis[i] = new SimpleContributionItem(knownMotifs[i]);
                }
                trainer.setSeedContributions(cg, seedCis);
            }
            trainer.setRestrictCrossover(restrictCrossover);
            trainer.setCrossOverProb(crossOverProb);
            trainer.setReplaceComponentProb(replaceComponentProb);
            needsInit = true;
        }
        
        EvaluationManager em;
        if (distributed) {
            DistributedEvaluationManager dem = new DistributedEvaluationManager(port, throughtputMonitor, crabDebug);
            dem.setLocalThreads(workerThreads);
            if (crabSignifier > 0) {
            	dem.setCrabSignifier(crabSignifier);
            }
            if (crabRate > 0) {
            	dem.setCrabRate(crabRate);
            }
            dem.setDisableBalancer(disableBalancer);
            em = dem;
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
        
        int lastLoggedCycle = trainer.getCycle();
        long lastLoggedTime = System.currentTimeMillis();
        while (true) {
            Trainer.WeightedModel twm = trainer.next();
            int c = trainer.getCycle();

            double evidenceLowerBound =  NativeMath.addLog2(trainer.getAccumulatedEvidence(), trainer.getPriorResidue() + twm.getLikelihood());
            double evidenceDif = evidenceLowerBound - trainer.getAccumulatedEvidence();
            
            if (sampleFile != null && (c % sampleInterval == 0)) {
                writeMotifs(trainer.getBestModel(), sampleFile + '.' + c + ".xms");
            }
            if ((c % logInterval) == 0) {
                long currentTime = System.currentTimeMillis();
	            System.out.print("" + c + '\t' + twm.getWeight() + '\t' + twm.getLikelihood());
	            if (cardinalityOut) {
		            for (int m = 0; m < numMotifs; ++m) {
		                System.out.print("\t" + countMix(twm.getModel(), m));
		            }
	            }
                System.out.print("\t" + trainer.getAccumulatedEvidence());
	            // System.out.print("\t" + trainer.getLikelihoodIQR());
	            System.out.print("\t" + (1.0 * (currentTime - lastLoggedTime)) / (c - lastLoggedCycle));
                // System.out.print("\t" + trainer.getAccumulatedEvidence());
                // System.out.print("\t" + trainer.getPriorResidue());
                // System.out.print("\t" + evidenceLowerBound);
                // System.out.print("\t" + (evidenceLowerBound - trainer.getAccumulatedEvidence()));
	            System.out.println();
	            lastLoggedTime = currentTime;
	            lastLoggedCycle = c;
            }
            
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
            
            if ((maxCycles > 0 && c >= maxCycles) || (terminateOnConvergance && evidenceDif < 0.01)) {
                if (outFile != null) {
                    writeMotifs(trainer.getBestModel(), outFile);
                }
                
                if (em instanceof DistributedEvaluationManager) {
                    ((DistributedEvaluationManager) em).shutdown();
                }
		System.exit(0); // hackity hack.
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
        	ContributionItem ci = model.getContribution(motifCG, m);
        	WeightMatrix wm = (WeightMatrix) ci.getItem();
        	
            Motif motif = new Motif();
            motif.setName("motif" + m);
            motif.setWeightMatrix(wm);
            // if (wm instanceof NMWeightMatrix) {
            //	motif.getAnnotation().setProperty("nmica.history_thread", "" + ((NMWeightMatrix) wm).trackingId());
            // }
            StringBuilder historyBuilder = new StringBuilder();
            Iterator<HistoryThread> li = ci.getHistoryThread().lineage();
            while (li.hasNext()) {
            	historyBuilder.append(li.next().getId());
            	if (li.hasNext()) {
            		historyBuilder.append(',');
            	}
            }
            motif.getAnnotation().setProperty("nmica.history_thread", historyBuilder.toString());
            motifList.add(motif);
        }
        Motif[] motifs = motifList.toArray(new Motif[0]);
        
        Annotation props = new SmallAnnotation();
        props.setProperty("creator.name", "nminfer");
        props.setProperty("creator.version", Version.VERSION);
        props.setProperty("input", seqDBfiles != null ? seqDBfiles[0].getName() : "*unknown*");
        props.setProperty("date", new SimpleDateFormat().format(new Date()));
        
        FileOutputStream fos = new FileOutputStream(fileName);
        MotifIOTools.writeMotifSetXML(fos, motifs, props);
        fos.close();
    }
    
    private int countMix(MultiICAModel m, int comp) {
        int c = 0;
        int d = m.getDataSet().length;
        for (int i = 0; i < d; ++i) {
            if (m.getMixture(i).get(comp) != 0.0) {
                ++c;
            }
        }
        return c;
    }
    
    private static SequenceDB loadDB(File f, FiniteAlphabet alphabet)
        throws Exception
    {
        SequenceIterator si = SeqIOTools.readFasta(new BufferedReader(new FileReader(f)),alphabet.getTokenization("token"));
        SequenceDB seqDB = new HashSequenceDB();
        while (si.hasNext()) {
            Sequence seq = si.nextSequence();
            if (seqDB.ids().contains(seq.getName())) {
                throw new IllegalIDException("Duplicate sequence name '" + seq.getName() + "'");
            }
            seqDB.addSequence(seq);
        }
        return seqDB;
    }
}
