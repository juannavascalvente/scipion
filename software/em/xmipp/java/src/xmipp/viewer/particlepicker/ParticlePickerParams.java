/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package xmipp.viewer.particlepicker;

import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.particlepicker.training.model.Mode;

/**
 *
 * @author airen
 */
public class ParticlePickerParams {
    
    private Options options;

    public final static String SCIPIONOPT = "scipion";
    public final static String INPUTOPT = "input";
    public final static String OUTPUTOPT = "output";
    public final static String MODEOPT = "mode";
    public final static String THREADSOPT = "threads";
    public final static String FASTOPT = "fast";
    public final static String INCOREOPT = "incore";
    private CommandLine cmdLine;
    
    
    public String inputfile;
    public String outputdir;
    public Integer port;
    public Mode mode;
    public Integer threads = 1;
    public boolean fast;
    public boolean incore;
    
    public ParticlePickerParams(Mode mode)
    {
    	this.mode = mode;
    }
    
    public ParticlePickerParams(String[] args)
    {
        try {
            defineArgs();
            processArgs(args);
        } catch (ParseException ex) {
            Logger.getLogger(ParticlePickerParams.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public void defineArgs() {
        options = new Options();
        options.addOption(INPUTOPT, true, "");
        options.addOption(OUTPUTOPT, true, "");
        options.addOption(MODEOPT, true, "");
        options.addOption(THREADSOPT, true, "");
        options.addOption(FASTOPT, true, "");
        options.addOption(INCOREOPT, true, "");
        options.addOption(SCIPIONOPT, true, "");
    }

    public void processArgs(String args[]) throws ParseException {

        String[] cmdargs;
        BasicParser parser = new BasicParser();
        cmdLine = parser.parse(options, args);
        inputfile = cmdLine.getOptionValue(INPUTOPT);
        outputdir = cmdLine.getOptionValue(OUTPUTOPT);
        mode = Mode.Manual;
        if(cmdLine.hasOption(MODEOPT))
        {
	        String str = cmdLine.getOptionValue(MODEOPT);
	        mode = Mode.getMode(str);
	        if(!(mode == Mode.Review || mode == Mode.ReadOnly))
	        	throw new IllegalArgumentException("Only Review and ReadOnly modes can be specified from the command line");
        }
        if (cmdLine.hasOption(THREADSOPT)) 
            threads = Integer.parseInt(cmdLine.getOptionValue(THREADSOPT));
        if (cmdLine.hasOption(FASTOPT)) 
            fast = Boolean.parseBoolean(cmdLine.getOptionValue(FASTOPT));
        if (cmdLine.hasOption(INCOREOPT)) 
            incore = Boolean.parseBoolean(cmdLine.getOptionValue(INCOREOPT));
        
       
        if (cmdLine.hasOption(SCIPIONOPT)) {
            cmdargs = cmdLine.getOptionValues(SCIPIONOPT);
            if(cmdargs != null)
                port = Integer.parseInt(cmdargs[0]);
        }

    }
    
    
    public boolean isScipion()
    {
        return cmdLine.hasOption(SCIPIONOPT);
    }
    
}
