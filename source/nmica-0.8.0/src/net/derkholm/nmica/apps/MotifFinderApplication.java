package net.derkholm.nmica.apps;

import net.derkholm.nmica.utils.CliTools;

public class MotifFinderApplication {
    public static void main(String[] args) throws Exception {
        MotifFinder app = new MotifFinder();
        args = CliTools.configureBean(app, args);
        app.main(args);
    }
}
