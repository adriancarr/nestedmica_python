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
 * Created on Aug 24, 2004
 */
package net.derkholm.nmica.demos;

import java.rmi.registry.LocateRegistry;
import java.rmi.registry.Registry;

import net.derkholm.nmica.utils.CliTools;

/**
 * @author thomas
 */
public class RMIObjectObjectPingClient {
    private String registryServer;
    private int registryPort = 1099;
    private String runKey = "pingServer";
    private int pings = 1000;

    public void setRegistryPort(int registryPort) {
        this.registryPort = registryPort;
    }
    public void setRegistryServer(String registryServer) {
        this.registryServer = registryServer;
    }
    public void setRunKey(String runKey) {
        this.runKey = runKey;
    }
    public void setPings(int p) {
        this.pings = p;
    }
    
    public static void main(String[] args) 
    		throws Exception
    {
        RMIObjectObjectPingClient app = new RMIObjectObjectPingClient();
        args = CliTools.configureBean(app, args);
        app.run(args);
    }
    
    public void run(String[] args)
    		throws Exception
    {
        Registry rmiReg  = LocateRegistry.getRegistry(registryServer, registryPort);
        Pingable responder = (Pingable) rmiReg.lookup(runKey);
        
        long before = System.currentTimeMillis();
        for (int p = 0; p < pings; ++p) {
            PingPacket pp = responder.objectObjectPing(new PingPacket(p, 5));
        }
        long after = System.currentTimeMillis();
        
        System.out.println("RMI roundtrip takes " + ((1.0 * (after - before)) / pings) + "ms");
    }
}
