within Modelica_Fluid;
package Interfaces 
  "Interfaces for steady state and unsteady, mixed-phase, multi-substance, incompressible and compressible flow" 
  
  annotation (Documentation(info="<html>
 
</html>", revisions="<html>
<ul>
<li><i>May 30, 2007</i>
       by Christoph Richter: moved everything back to its original position in Modelica_Fluid.</li>
<li><i>Apr. 20, 2007</i>
       by Christoph Richter: moved parts of the original package from Modelica_Fluid
       to the development branch of Modelica 2.2.2.</li>
<li><i>Nov. 2, 2005</i>
       by Francesco Casella: restructured after 45th Design Meeting.</li>
<li><i>Nov. 20-21, 2002</i>
       by Hilding Elmqvist, Mike Tiller, Allan Watson, John Batteh, Chuck Newman,
       Jonas Eborn: Improved at the 32nd Modelica Design Meeting.
<li><i>Nov. 11, 2002</i>
       by Hilding Elmqvist, Martin Otter: improved version.</li>
<li><i>Nov. 6, 2002</i>
       by Hilding Elmqvist: first version.</li>
<li><i>Aug. 11, 2002</i>
       by Martin Otter: Improved according to discussion with Hilding
       Elmqvist and Hubertus Tummescheit.<br>
       The PortVicinity model is manually
       expanded in the base models.<br>
       The Volume used for components is renamed
       PartialComponentVolume.<br>
       A new volume model \"Fluid.Components.PortVolume\"
       introduced that has the medium properties of the port to which it is
       connected.<br>
       Fluid.Interfaces.PartialTwoPortTransport is a component
       for elementary two port transport elements, whereas PartialTwoPort
       is a component for a container component.</li>
</li>
</ul>
</html>"));
  
  extends Modelica.Icons.Library;
  
  connector FluidPort 
    "Interface for quasi one-dimensional fluid flow in a piping network (incompressible or compressible, one or more phases, one or more substances)" 
    
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
      "Medium model" annotation (choicesAllMatching=true);
    
    flow Medium.MassFlowRate m_flow 
      "Mass flow rate from the connection point into the component";
    Medium.AbsolutePressure p "Pressure in the connection point";
    Medium.SpecificEnthalpy h_outflow 
      "Specific enthalpy close to the connection point if m_flow < 0" 
      annotation(stream);
    Medium.MassFraction Xi_outflow[Medium.nXi] 
      "Independent mixture mass fractions m_i/m close to the connection point if m_flow < 0"
      annotation(stream);
    Medium.ExtraProperty C_outflow[Medium.nC] 
      "Properties c_i/m close to the connection point if m_flow < 0" 
      annotation(stream);
  end FluidPort;
  
  connector FluidPort_a "Generic fluid connector at design inlet" 
    extends FluidPort;
    annotation (defaultComponentName="port_a",
                Diagram(Ellipse(extent=[-40,40; 40,-40], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=69,
            rgbfillColor={0,127,255})),
                               Text(extent=[-150,110; 150,50],   string="%name")),
         Icon(Ellipse(extent=[-100, 100; 100, -100], style(color=69,
              fillColor=69)), Ellipse(extent=[-100, 100; 100, -100], style(color=16,
              fillColor=69))));
  end FluidPort_a;
  
  connector FluidPort_b "Generic fluid connector at design outlet" 
    extends FluidPort;
    annotation (defaultComponentName="port_b",
                Diagram(Ellipse(extent=[-40,40; 40,-40], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=69,
            rgbfillColor={0,127,255})),
                               Ellipse(extent=[-30,30; 30,-30],   style(color=69,
               fillColor=7)), Text(extent=[-150,110; 150,50],   string="%name")),
         Icon(Ellipse(extent=[-100, 100; 100, -100], style(color=69,
              fillColor=69)), Ellipse(extent=[-100, 100; 100, -100], style(color=16,
              fillColor=69)), Ellipse(extent=[-80, 80; 80, -80], style(color=69,
               fillColor=7))));
  end FluidPort_b;
  
  connector FluidStatePort_a 
    "Fluid connector at design inlet with potential pressure state" 
    extends FluidPort;
    annotation (defaultComponentName="port_a",
                Diagram(Ellipse(extent=[-40,40; 40,-40], style(
            color=3,
            fillColor=69,
            rgbfillColor={0,127,255})),
                               Text(extent=[-150,110; 150,50],   string="%name")),
         Icon(Ellipse(extent=[-100, 100; 100, -100], style(color=69,
              fillColor=69)), Ellipse(extent=[-100, 100; 100, -100], style(color=16,
              fillColor=69)),
              Ellipse(extent=[-18,20; 22,-20], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0}))));
  end FluidStatePort_a;
  
  connector FluidStatePort_b 
    "Fluid connector at design outlet with potential pressure state" 
   extends FluidPort;
    annotation (defaultComponentName="port_b",
                Diagram(Ellipse(extent=[-40,40; 40,-40], style(
            color=3,
            fillColor=69,
            rgbfillColor={0,127,255})),
                               Ellipse(extent=[-30,30; 30,-30],   style(color=69,
               fillColor=7)), Text(extent=[-150,110; 150,50],   string="%name")),
         Icon(Ellipse(extent=[-100, 100; 100, -100], style(color=69,
              fillColor=69)), Ellipse(extent=[-100, 100; 100, -100], style(color=16,
              fillColor=69)), Ellipse(extent=[-80, 80; 80, -80], style(color=69,
               fillColor=7)),
              Ellipse(extent=[-18,20; 22,-20], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0}))));
  end FluidStatePort_b;
  
  connector FluidPorts_a 
    "Fluid connector with filled, large icon to be used for vectors of FluidPorts (vector dimensions must be added after dragging)" 
    extends FluidPort;
    annotation (defaultComponentName="ports_a",
                Diagram(       Text(extent=[-75,130; 75,100],  string="%name"),
        Rectangle(extent=[-25,100; 25,-100], style(
            color=69,
            fillColor=7,
            rgbfillColor={255,255,255})),
        Ellipse(extent=[-25,90; 25,40],style(color=16,fillColor=69)),
        Ellipse(extent=[-25,25; 25,-25],style(color=16,fillColor=69)),
        Ellipse(extent=[-25,-40; 25,-90], style(color=16,fillColor=69))),
         Icon(
        Rectangle(extent=[-50,200; 50,-200], style(
            color=69,
            fillColor=7,
            rgbfillColor={255,255,255})),
                              Ellipse(extent=[-50,180; 50,80],       style(color=16,
              fillColor=69)), Ellipse(extent=[-50,50; 50,-50],       style(color=16,
              fillColor=69)), Ellipse(extent=[-50,-80; 50,-180],     style(color=16,
              fillColor=69))),
      Coordsys(
        extent=[-50,-200; 50,200],
        grid=[1,1],
        scale=0.2));
  end FluidPorts_a;
  
  connector FluidPorts_b 
    "Fluid connector with outlined, large icon to be used for vectors of FluidPorts (vector dimensions must be added after dragging)" 
    extends FluidPort;
    annotation (defaultComponentName="ports_b",
                Diagram(       Text(extent=[-75,130; 75,100],  string="%name"),
        Rectangle(extent=[-25,100; 25,-100], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255})),
        Ellipse(extent=[-25,90; 25,40],style(color=16,fillColor=69)),
        Ellipse(extent=[-25,25; 25,-25],style(color=16,fillColor=69)),
        Ellipse(extent=[-25,-40; 25,-90], style(color=16,fillColor=69)),
        Ellipse(extent=[-15,-50; 15,-80], style(color=69, fillColor=7)),
        Ellipse(extent=[-15,15; 15,-15], style(color=69, fillColor=7)),
        Ellipse(extent=[-15,50; 15,80], style(color=69, fillColor=7))),
         Icon(
        Rectangle(extent=[-50,200; 50,-200], style(
            color=69,
            fillColor=7,
            rgbfillColor={255,255,255})),
                              Ellipse(extent=[-50,180; 50,80],       style(color=16,
              fillColor=69)), Ellipse(extent=[-50,50; 50,-50],       style(color=16,
              fillColor=69)), Ellipse(extent=[-50,-80; 50,-180],     style(color=16,
              fillColor=69)),
            Ellipse(extent=[-30,30; 30,-30], style(color=69, fillColor=7)),
            Ellipse(extent=[-30,100; 30,160], style(color=69, fillColor=7)),
            Ellipse(extent=[-30,-100; 30,-160], style(color=69, fillColor=7))),
      Coordsys(
        extent=[-50,-200; 50,200],
        grid=[1,1],
        scale=0.2));
  end FluidPorts_b;
  
  connector FluidStatePorts_a 
    "Fluid connector with filled, large icon to be used for vectors of FluidPorts (vector dimensions must be added after dragging)" 
    extends FluidPort;
    annotation (defaultComponentName="ports_a",
                Diagram(       Text(extent=[-75,130; 75,100],  string="%name"),
        Rectangle(extent=[-25,100; 25,-100], style(
            color=69,
            fillColor=7,
            rgbfillColor={255,255,255})),
        Ellipse(extent=[-25,90; 25,40],style(color=16,fillColor=69)),
        Ellipse(extent=[-25,25; 25,-25],style(color=16,fillColor=69)),
        Ellipse(extent=[-25,-40; 25,-90], style(color=16,fillColor=69))),
         Icon(
        Rectangle(extent=[-50,200; 50,-200], style(
            color=69,
            fillColor=7,
            rgbfillColor={255,255,255})),
                              Ellipse(extent=[-50,180; 50,80],       style(color=16,
              fillColor=69)), Ellipse(extent=[-50,50; 50,-50],       style(color=16,
              fillColor=69)), Ellipse(extent=[-50,-80; 50,-180],     style(color=16,
              fillColor=69)),
        Ellipse(extent=[-20,150; 20,110], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0})),
        Ellipse(extent=[-20,20; 20,-20], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0})),
        Ellipse(extent=[-19,-111; 21,-151], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0}))),
      Coordsys(
        extent=[-50,-200; 50,200],
        grid=[1,1],
        scale=0.2));
  end FluidStatePorts_a;
  
  connector FluidStatePorts_b 
    "Fluid connector with outlined, large icon to be used for vectors of FluidPorts (vector dimensions must be added after dragging)" 
    extends FluidPort;
    annotation (defaultComponentName="ports_b",
                Diagram(       Text(extent=[-75,130; 75,100],  string="%name"),
        Rectangle(extent=[-25,100; 25,-100], style(
            color=69,
            fillColor=7,
            rgbfillColor={255,255,255})),
        Ellipse(extent=[-25,90; 25,40],style(color=16,fillColor=69)),
        Ellipse(extent=[-25,25; 25,-25],style(color=16,fillColor=69)),
        Ellipse(extent=[-25,-40; 25,-90], style(color=16,fillColor=69)),
        Ellipse(extent=[-15,-50; 15,-80], style(color=69, fillColor=7)),
        Ellipse(extent=[-15,15; 15,-15], style(color=69, fillColor=7)),
        Ellipse(extent=[-15,50; 15,80], style(color=69, fillColor=7))),
         Icon(
        Rectangle(extent=[-50,200; 50,-200], style(
            color=69,
            fillColor=7,
            rgbfillColor={255,255,255})),
                              Ellipse(extent=[-50,180; 50,80],       style(color=16,
              fillColor=69)), Ellipse(extent=[-50,50; 50,-50],       style(color=16,
              fillColor=69)), Ellipse(extent=[-50,-80; 50,-180],     style(color=16,
              fillColor=69)),
            Ellipse(extent=[-45,44; 45,-44], style(color=69, fillColor=7)),
            Ellipse(extent=[-43,88; 43,172],  style(color=69, fillColor=7)),
            Ellipse(extent=[-43,-84; 45,-175],  style(color=69, fillColor=7)),
        Ellipse(extent=[-20,151; 20,111], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0})),
        Ellipse(extent=[-20,21; 20,-19], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0})),
        Ellipse(extent=[-19,-110; 21,-150], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0}))),
      Coordsys(
        extent=[-50,-200; 50,200],
        grid=[1,1],
        scale=0.2));
  end FluidStatePorts_b;
  
  connector HeatPorts_a 
    "HeatPort connector with filled, large icon to be used for vectors of HeatPorts (vector dimensions must be added after dragging)" 
    extends Modelica.Thermal.HeatTransfer.Interfaces.HeatPort;
    annotation (defaultComponentName="heatPorts_a",
                Diagram(       Text(extent=[-75,130; 75,100],  string="%name"),
        Rectangle(extent=[-25,100; 25,-100], style(
            color=42,
            rgbcolor={127,0,0},
            fillColor=7,
            rgbfillColor={255,255,255})),
        Rectangle(extent=[-21,89; 22,44], style(
            color=42,
            rgbcolor={127,0,0},
            fillColor=42,
            rgbfillColor={127,0,0},
            fillPattern=1)),
        Rectangle(extent=[-21,22; 22,-23], style(
            color=42,
            rgbcolor={127,0,0},
            fillColor=42,
            rgbfillColor={127,0,0},
            fillPattern=1)),
        Rectangle(extent=[-21.5,-43; 21.5,-88], style(
            color=42,
            rgbcolor={127,0,0},
            fillColor=42,
            rgbfillColor={127,0,0},
            fillPattern=1))),
         Icon(
        Rectangle(extent=[-50,200; 50,-200], style(
            color=42,
            rgbcolor={127,0,0},
            fillColor=7,
            rgbfillColor={255,255,255})),
        Rectangle(extent=[-44,176; 44,86], style(
            color=42,
            rgbcolor={127,0,0},
            fillColor=42,
            rgbfillColor={127,0,0},
            fillPattern=1)),
        Rectangle(extent=[-43,46; 45,-44], style(
            color=42,
            rgbcolor={127,0,0},
            fillColor=42,
            rgbfillColor={127,0,0},
            fillPattern=1)),
        Rectangle(extent=[-44,-86; 44,-176], style(
            color=42,
            rgbcolor={127,0,0},
            fillColor=42,
            rgbfillColor={127,0,0},
            fillPattern=1))),
      Coordsys(
        extent=[-50,-200; 50,200],
        grid=[1,1],
        scale=0.2));
  end HeatPorts_a;
  
  connector HeatPorts_b 
    "HeatPort connector with filled, large icon to be used for vectors of HeatPorts (vector dimensions must be added after dragging)" 
    extends Modelica.Thermal.HeatTransfer.Interfaces.HeatPort;
    annotation (defaultComponentName="heatPorts_b",
                Diagram(       Text(extent=[-75,130; 75,100],  string="%name"),
        Rectangle(extent=[-25,100; 25,-100], style(
            color=42,
            rgbcolor={127,0,0},
            fillColor=7,
            rgbfillColor={255,255,255})),
        Rectangle(extent=[-21,88; 22,43], style(
            color=42,
            rgbcolor={127,0,0},
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Rectangle(extent=[-21,22; 22,-23], style(
            color=42,
            rgbcolor={127,0,0},
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Rectangle(extent=[-21,-43; 22,-88], style(
            color=42,
            rgbcolor={127,0,0},
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1))),
         Icon(
        Rectangle(extent=[-50,200; 50,-200], style(
            color=42,
            rgbcolor={127,0,0},
            fillColor=7,
            rgbfillColor={255,255,255})),
        Rectangle(extent=[-43,175; 45,85], style(
            color=42,
            rgbcolor={127,0,0},
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Rectangle(extent=[-44,46; 44,-44], style(
            color=42,
            rgbcolor={127,0,0},
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Rectangle(extent=[-44,-86; 44,-176], style(
            color=42,
            rgbcolor={127,0,0},
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1))),
      Coordsys(
        extent=[-50,-200; 50,200],
        grid=[1,1],
        scale=0.2));
  end HeatPorts_b;
end Interfaces;
