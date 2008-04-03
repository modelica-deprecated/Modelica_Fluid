within FluidSandbox;
package Interfaces 
  "Interfaces to implement different elements of a numerical solution method for fluid dynamics" 
  
  partial package PartialFluidInterface 
    "Partial fluid interface implementation" 
    extends Icons.FluidInterface;
    
    constant Boolean usesNewConnectionSemantics=false 
      "Does this approach use yet unsupported connection semantics?";
    
    replaceable package FluidDiscretization = PartialFluidDiscretization 
      "Fluid discretization currently used - useful to adapt interfaces to discretization e.g. momentum flow";
    
    replaceable partial connector FluidPort 
      "Partial connector for FluidInterface" 
      
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Medium model" annotation (choicesAllMatching=true);
      
    end FluidPort;
    
    replaceable partial connector FluidPort_a 
      "Generic fluid connector at design inlet" 
      extends FluidPort;
      annotation (
        defaultComponentName="port_a",
        Diagram(Ellipse(extent=[-40,40; 40,-40], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=69,
              rgbfillColor={0,127,255})), Text(extent=[-150,110; 150,50],
              string="%name")),
        Icon(Ellipse(extent=[-100,100; 100,-100], style(color=69, fillColor=69)),
            Ellipse(extent=[-100,100; 100,-100], style(color=16, fillColor=69))));
    end FluidPort_a;
    
    replaceable partial connector FluidPort_b 
      "Generic fluid connector at design outlet" 
      extends FluidPort;
      annotation (
        defaultComponentName="port_b",
        Diagram(
          Ellipse(extent=[-40,40; 40,-40], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=69,
              rgbfillColor={0,127,255})),
          Ellipse(extent=[-30,30; 30,-30], style(color=69, fillColor=7)),
          Text(extent=[-150,110; 150,50], string="%name")),
        Icon(
          Ellipse(extent=[-100,100; 100,-100], style(color=69, fillColor=69)),
          Ellipse(extent=[-100,100; 100,-100], style(color=16, fillColor=69)),
          Ellipse(extent=[-80,80; 80,-80], style(color=69, fillColor=7))));
    end FluidPort_b;
    
    replaceable partial model ConnectionSemantics 
      "Allows to implement unsupported connection semantics" 
      
      // Medium model
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Medium model" annotation (choicesAllMatching=true);
      
      // Interfaces
      FluidPort_a port_a(redeclare package Medium = Medium) "Fluid inlet port" 
                           annotation (extent=[-110,-10; -90,10]);
      FluidPort_b port_b(redeclare package Medium = Medium) "Fluid outlet port"
                            annotation (extent=[90,-10; 110,10]);
      
      annotation (
        defaultComponentName="semantics",
        Coordsys(extent=[-100,-20; 100,20], scale=0.05),
        Icon(Text(extent=[-3,0; 3,-64], string=" ")),
        Diagram);
    end ConnectionSemantics;
    
    replaceable partial model PartialLumpedVolume 
      "Mixing volume without connectors" 
      
    // Medium model
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Medium in the component" 
          annotation (choicesAllMatching = true);
      
      // Interfaces
      FluidPort_a port_a[n_a](redeclare package Medium = Medium) 
        "Fluid inlet port" annotation (extent=[-110,-10; -90,10]);
      FluidPort_a port_b[n_b](redeclare package Medium = Medium) 
        "Fluid outlet port (Port_a, too!)" annotation (extent=[90,-10; 110,10]);
      parameter Integer n_a=1 "Number of port_a's";
      parameter Integer n_b=1 "Number of port_b's";
      
    // BaseProperties instance
      Medium.BaseProperties medium 
        "No start or guess values as states have to be chosen in submodel";
      
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort 
        annotation (extent=[-10,-110; 10,-90]);
      
      // add a modifier to these declaration equations to override this
      Medium.SpecificEnthalpy h_a_outflow = medium.h 
        "Specific enthalpy at port_a if fluid exits";
      Medium.MassFraction Xi_a_outflow[Medium.nXi] = medium.Xi 
        "Independent mixture mass fractions m_i/m in the connection point a if fluid exits";
      
      Medium.SpecificEnthalpy h_b_outflow = medium.h 
        "Specific enthalpy at port_a if fluid exits";
      Medium.MassFraction Xi_b_outflow[Medium.nXi] = medium.Xi 
        "Independent mixture mass fractions m_i/m in the connection point b if fluid exits";
      
      Medium.MassFlowRate m_flow_net "Net mass flow into the volume";
      Medium.MassFlowRate mXi_flow_net[Medium.nXi] 
        "Net substance mass flow into the volume";
      Medium.EnthalpyFlowRate H_flow_net "Net enthalpy flow into the volume";
      
      // sensors /////////////////  
      parameter Boolean provide_p = false "Provide pressure?" annotation(Evaluate=true, Dialog(descriptionLabel=true, tab="Sensors"));
      parameter Boolean provide_T = false "Provide temperature?" annotation(Evaluate=true, Dialog(descriptionLabel=true, tab="Sensors"));
      
      Modelica.Blocks.Interfaces.RealOutput p_sensor(redeclare type SignalType 
          = SI.Pressure) if provide_p annotation (extent=[-120,70; -100,90], rotation=180);
      Modelica.Blocks.Interfaces.RealOutput T_sensor(redeclare type SignalType 
          = SI.Temperature) if provide_T annotation (extent=[-120,40; -100,60], rotation=180);
    protected 
      Modelica.Blocks.Interfaces.RealOutput calc_p(redeclare type SignalType = 
            SI.Pressure) annotation (extent=[-74,40; -66,48], rotation=90);
      Modelica.Blocks.Interfaces.RealOutput calc_T(redeclare type SignalType = 
            SI.Temperature) annotation (extent=[-84,40; -76,48], rotation=90);
    equation 
      // boundary conditions heat port
      heatPort.T = medium.T;
      
      // sensors
      calc_p = medium.p;
      calc_T = medium.T;
      connect(T_sensor,calc_T)  annotation (points=[-110,50; -80,50; -80,44], style(
          color=74,
          rgbcolor={0,0,127},
          smooth=0));
      connect(calc_p,p_sensor)  annotation (points=[-70,44; -70,80; -110,80], style(
          color=74,
          rgbcolor={0,0,127},
          smooth=0));
      
      annotation (
        Icon(Text(extent=[-144,178; 146,116], string="%name"), Text(
            extent=[-130,-108; 144,-150],
            style(color=0),
            string="V=%V")),
        Documentation(info="<html>
Base class for an ideally mixed fluid volume with two ports and the ability to store mass and energy. The following source terms are part of the energy balance and must be specified in the extending class:
<ul>
<li><tt>Qs_flow</tt>, e.g. convective or latent heat flow rate across segment boundary, and</li> <li><tt>Ws_flow</tt>, work term, e.g. p*der(V) if the volume is not constant</li>
</ul>
The component volume <tt>V_lumped</tt> is also a variable which needs to be set in the extending class to complete the model.
</html>"),
        Diagram);
      
    end PartialLumpedVolume;
    
    replaceable partial model PartialTransport 
      "Partial isenthalpic element transporting fluid between two ports without storing mass or energy (two Port_b's)" 
      
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Medium in the component"                                          annotation (
          choicesAllMatching =                                                                            true);
      
      Medium.MassFlowRate m_flow 
        "Mass flow rate from port_a to port_b (m_flow > 0 is design flow direction)";
      SI.Pressure dp "Pressure drop between port_a and port_b";
      
      // Used for some approaches only
      parameter SI.ThermalConductance G=0 
        "Conductive heat flow through component (used in some fluid interfaces only)";
      parameter SI.DiffusionCoefficient H=0 
        "Mass diffusion through component (used in some fluid interfaces only)";
      
      // Properties for fluid flowing in design direction, have to be provided by PartialTwoPortTransport
      // Used in pressure drop correlation
      Medium.AbsolutePressure p_designDirection 
        "Ideally: Upstream p if fluid flew in design direction (even if currently not the case). At least: Upstream p if fluid actually flows in design direction.";
      Medium.SpecificEnthalpy h_designDirection 
        "Ideally: Upstream h if fluid flew in design direction (even if currently not the case). At least: Upstream h if fluid actually flows in design direction.";
      Medium.MassFraction Xi_designDirection[Medium.nXi] 
        "Ideally: Upstream Xi if fluid flew in design direction (even if currently not the case). At least: Upstream Xi if fluid actually flows in design direction.";
      
      // Properties for fluid flowing in non-design direction, have to be provided by PartialTwoPortTransport
      // Used in pressure drop correlation
      Medium.AbsolutePressure p_nonDesignDirection 
        "Ideally: Upstream p if fluid flew against design direction (even if currently not the case). At least: Upstream p if fluid actually flows against design direction.";
      Medium.SpecificEnthalpy h_nonDesignDirection 
        "Ideally: Upstream h if fluid flew against design direction (even if currently not the case). At least: Upstream h if fluid actually flows against design direction.";
      Medium.MassFraction Xi_nonDesignDirection[Medium.nXi] 
        "Ideally: Upstream Xi if fluid flew against design direction (even if currently not the case). At least: Upstream Xi if fluid actually flows against design direction.";
      
      // Sensor outputs //////////////////////////////////////////////////////////////////////////////
      Modelica.Blocks.Interfaces.RealOutput p_a(redeclare type SignalType = 
            SI.Pressure) if provide_p_a annotation (extent=[-120,70; -100,90], rotation=180);
      Modelica.Blocks.Interfaces.RealOutput p_b(redeclare type SignalType = 
            SI.Pressure) if provide_p_b annotation (extent=[100,70; 120,90], rotation=0);
      Modelica.Blocks.Interfaces.RealOutput T_a(redeclare type SignalType = 
            SI.Temperature) if provide_T_a annotation (extent=[-120,40; -100,60], rotation=180);
      Modelica.Blocks.Interfaces.RealOutput T_b(redeclare type SignalType = 
            SI.Temperature) if provide_T_b annotation (extent=[100,40; 120,60], rotation=0);
      Modelica.Blocks.Interfaces.RealOutput m_flow_ab(redeclare type SignalType
          = SI.MassFlowRate) if provide_m_flow_ab annotation (extent=[-70,100; -50,
            120],                                                                       rotation=90);
      
      parameter Boolean provide_p_a = false "Provide pressure at port a?" annotation(Evaluate=true, Dialog(descriptionLabel=true, tab="Sensors"));
      parameter Boolean provide_p_b = false "Provide pressure at port b?" annotation(Evaluate=true, Dialog(descriptionLabel=true, tab="Sensors"));
      parameter Boolean provide_T_a = false "Provide temperature at port a?" annotation(Evaluate=true, Dialog(descriptionLabel=true, tab="Sensors"));
      parameter Boolean provide_T_b = false "Provide temperature at port b?" annotation(Evaluate=true, Dialog(descriptionLabel=true, tab="Sensors"));
      parameter Boolean provide_m_flow_ab = false "Provide mass flow rate?" annotation(Evaluate=true, Dialog(descriptionLabel=true, tab="Sensors"));
      
    protected 
      Modelica.Blocks.Interfaces.RealOutput calc_T_a(redeclare type SignalType 
          = SI.Temperature) annotation (extent=[-84,40; -76,48], rotation=90);
      Modelica.Blocks.Interfaces.RealOutput calc_p_a(redeclare type SignalType 
          = SI.Pressure) annotation (extent=[-74,40; -66,48], rotation=90);
      Modelica.Blocks.Interfaces.RealOutput calc_m_flow_ab(redeclare type 
          SignalType = SI.MassFlowRate) 
        annotation (extent=[-64,40; -56,48], rotation=90);
      Modelica.Blocks.Interfaces.RealOutput calc_p_b(redeclare type SignalType 
          = SI.Pressure) annotation (extent=[-54,40; -46,48], rotation=90);
      Modelica.Blocks.Interfaces.RealOutput calc_T_b(redeclare type SignalType 
          = SI.Temperature) annotation (extent=[-44,40; -36,48], rotation=90);
      
      // Using Medium.temperature(Medium.setState_phX()) for temperature sensor results in numeric Jacobian; using BaseProperties instead
      Medium.BaseProperties calc_T_a_medium;
      Medium.BaseProperties calc_T_b_medium;
    equation 
      connect(T_a, calc_T_a) annotation (points=[-110,50; -80,50; -80,44], style(
          color=74,
          rgbcolor={0,0,127},
          smooth=0));
      connect(calc_p_a, p_a) annotation (points=[-70,44; -70,80; -110,80], style(
          color=74,
          rgbcolor={0,0,127},
          smooth=0));
      connect(m_flow_ab, calc_m_flow_ab) annotation (points=[-60,110; -60,44],
          style(
          color=74,
          rgbcolor={0,0,127},
          smooth=0));
      connect(calc_p_b, p_b) annotation (points=[-50,44; -50,80; 110,80], style(
          color=74,
          rgbcolor={0,0,127},
          smooth=0));
      connect(calc_T_b, T_b) annotation (points=[-40,44; -40,50; 110,50], style(
          color=74,
          rgbcolor={0,0,127},
          smooth=0));
      
      annotation (Icon(Text(
            extent=[-144,-61; 144,-110],
            string="%name",
            style(gradient=2, fillColor=69)),
          Text(
            extent=[-132,116; -100,100],
            style(color=3, rgbcolor={0,0,255}),
            string="p_a"),
          Text(
            extent=[100,116; 132,100],
            style(color=3, rgbcolor={0,0,255}),
            string="p_b"),
          Text(
            extent=[100,30; 132,14],
            style(color=3, rgbcolor={0,0,255}),
            string="T_b"),
          Text(
            extent=[-132,30; -100,14],
            style(color=3, rgbcolor={0,0,255}),
            string="T_a")), Diagram);
    end PartialTransport;
    
    replaceable partial model PartialTransportIsenthalpic 
      "Partial isenthalpic element transporting fluid between two ports without storing mass or energy (two Port_b's)" 
      
      extends PartialTransport;
      
      FluidPort_b port_a(redeclare package Medium = Medium) 
        "Fluid connector a (positive design flow direction is from port_a to port_b)"
        annotation (extent=[-110,-10; -90,10]);
      FluidPort_b port_b(redeclare package Medium = Medium) 
        "Fluid connector b (positive design flow direction is from port_a to port_b)"
        annotation (extent=[110,-10; 90,10]);
      
    end PartialTransportIsenthalpic;
    
    replaceable partial model PartialTransportIsenthalpicAA 
      "Partial isenthalpic element transporting fluid between two ports without storing mass or energy (two Port_a's, not supported for all interfaces)" 
      
      extends PartialTransport;
      
      FluidPort_a port_a(redeclare package Medium = Medium) 
        "Fluid connector a (positive design flow direction is from port_a to port_b)"
        annotation (extent=[-110,-10; -90,10]);
      FluidPort_a port_b(redeclare package Medium = Medium) 
        "Fluid connector b (positive design flow direction is from port_a to port_b)"
        annotation (extent=[110,-10; 90,10]);
      
    end PartialTransportIsenthalpicAA;
    
    replaceable partial model PartialTransportIsenthalpicAB 
      "Partial isenthalpic element transporting fluid between two ports without storing mass or energy (a Port_a and Port_b each, not supported for all interfaces)" 
      
      extends PartialTransport;
      
      FluidPort_a port_a(redeclare package Medium = Medium) 
        "Fluid connector a (positive design flow direction is from port_a to port_b)"
        annotation (extent=[-110,-10; -90,10]);
      FluidPort_b port_b(redeclare package Medium = Medium) 
        "Fluid connector b (positive design flow direction is from port_a to port_b)"
        annotation (extent=[110,-10; 90,10]);
      
    end PartialTransportIsenthalpicAB;
    
    replaceable partial model PartialTransportIsentropic 
      "Partial Isentropic element transporting fluid between two ports without storing mass or energy (two Port_b's)" 
      
      import Modelica_Fluid.Types;
      
      extends PartialTransport;
      
      parameter Types.FlowDirection.Temp flowDirection=Types.FlowDirection.Bidirectional 
        "Unidirectional (port_a -> port_b) or bidirectional flow component" 
         annotation(Evaluate=true, Dialog(tab="Advanced"));
      
      // Limit the decomposition complexity and repeat declarations
      parameter Real eta_ise = 1 "Isentropic efficiency";
      Modelica.SIunits.Power P_mechanical 
        "Power exchanged through work done on the environment";
      
      FluidPort_b port_a(redeclare package Medium = Medium) 
        "Fluid connector a (positive design flow direction is from port_a to port_b)"
          annotation (extent=[-110,-10; -90,10]);
      FluidPort_b port_b(redeclare package Medium = Medium) 
        "Fluid connector b (positive design flow direction is from port_a to port_b)"
          annotation (extent=[110,-10; 90,10]);
      
    protected 
      Medium.BaseProperties medium_designDirection 
        "Ideally: Upstream properties if fluid flew in design direction (even if currently not the case). At least: Upstream properties if fluid actually flows in design direction.";
      Medium.BaseProperties medium_nonDesignDirection 
        "Ideally: Upstream properties if fluid flew against design direction (even if currently not the case). At least: Upstream properties if fluid actually flows against design direction.";
      
    equation 
      // Media instances may only be used if sensible
      medium_designDirection.p = p_designDirection;
      medium_designDirection.h = h_designDirection;
      medium_designDirection.Xi = Xi_designDirection;
      medium_nonDesignDirection.p = if flowDirection==Types.FlowDirection.Bidirectional then p_nonDesignDirection else Medium.p_default;
      medium_nonDesignDirection.h = if flowDirection==Types.FlowDirection.Bidirectional then h_nonDesignDirection else Medium.h_default;
      medium_nonDesignDirection.Xi = if flowDirection==Types.FlowDirection.Bidirectional then Xi_nonDesignDirection else zeros(Medium.nXi);
      
      assert(flowDirection==Types.FlowDirection.Bidirectional or m_flow>-0.001, "PartialTransportIsentropic: Mass flow direction was said to be positive only but is not.");
      
    end PartialTransportIsentropic;
    
    replaceable partial model PartialTransportIsentropicAA 
      "Partial Isentropic element transporting fluid between two ports without storing mass or energy (two Port_a's, not supported for all interfaces)" 
      
      import Modelica_Fluid.Types;
      
      extends PartialTransport;
      
      parameter Types.FlowDirection.Temp flowDirection=Types.FlowDirection.Bidirectional 
        "Unidirectional (port_a -> port_b) or bidirectional flow component" 
         annotation(Evauate=true, Dialog(tab="Advanced"));
      
      // Limit the decomposition complexity and repeat declarations
      parameter Real eta_ise = 1 "Isentropic efficiency";
      Modelica.SIunits.Power P_mechanical 
        "Power exchanged through work done on the environment";
      
      FluidPort_a port_a(redeclare package Medium = Medium) 
        "Fluid connector a (positive design flow direction is from port_a to port_b)"
        annotation (extent=[-110,-10; -90,10]);
      FluidPort_a port_b(redeclare package Medium = Medium) 
        "Fluid connector b (positive design flow direction is from port_a to port_b)"
        annotation (extent=[110,-10; 90,10]);
      
    protected 
      Medium.BaseProperties medium_designDirection 
        "Ideally: Upstream properties if fluid flew in design direction (even if currently not the case). At least: Upstream properties if fluid actually flows in design direction.";
      Medium.BaseProperties medium_nonDesignDirection 
        "Ideally: Upstream properties if fluid flew against design direction (even if currently not the case). At least: Upstream properties if fluid actually flows against design direction.";
      
    equation 
      // Media instances may only be used if sensible
      medium_designDirection.p = p_designDirection;
      medium_designDirection.h = h_designDirection;
      medium_designDirection.Xi = Xi_designDirection;
      medium_nonDesignDirection.p = if flowDirection==Types.FlowDirection.Bidirectional then p_nonDesignDirection else Medium.p_default;
      medium_nonDesignDirection.h = if flowDirection==Types.FlowDirection.Bidirectional then h_nonDesignDirection else Medium.h_default;
      medium_nonDesignDirection.Xi = if flowDirection==Types.FlowDirection.Bidirectional then Xi_nonDesignDirection else zeros(Medium.nXi);
      
      assert(flowDirection==Types.FlowDirection.Bidirectional or m_flow>-0.001, "PartialTransportIsentropicAA: Mass flow direction was said to be positive only but is not.");
      
    end PartialTransportIsentropicAA;
    
    replaceable partial model PartialTransportIsentropicAB 
      "Partial Isentropic element transporting fluid between two ports without storing mass or energy (a Port_a and Port_b each, not supported for all interfaces)" 
      
      import Modelica_Fluid.Types;
      
      extends PartialTransport;
      
      parameter Types.FlowDirection.Temp flowDirection=Types.FlowDirection.Bidirectional 
        "Unidirectional (port_a -> port_b) or bidirectional flow component" 
         annotation(Evauate=true, Dialog(tab="Advanced"));
      
      // Limit the decomposition complexity and repeat declarations
      parameter Real eta_ise = 1 "Isentropic efficiency";
      Modelica.SIunits.Power P_mechanical 
        "Power exchanged through work done on the environment";
      
      FluidPort_a port_a(redeclare package Medium = Medium) 
        "Fluid connector a (positive design flow direction is from port_a to port_b)"
        annotation (extent=[-110,-10; -90,10]);
      FluidPort_b port_b(redeclare package Medium = Medium) 
        "Fluid connector b (positive design flow direction is from port_a to port_b)"
        annotation (extent=[110,-10; 90,10]);
      
    protected 
      Medium.BaseProperties medium_designDirection 
        "Ideally: Upstream properties if fluid flew in design direction (even if currently not the case). At least: Upstream properties if fluid actually flows in design direction.";
      Medium.BaseProperties medium_nonDesignDirection 
        "Ideally: Upstream properties if fluid flew against design direction (even if currently not the case). At least: Upstream properties if fluid actually flows against design direction.";
      
    equation 
      // Media instances may only be used if sensible
      medium_designDirection.p = p_designDirection;
      medium_designDirection.h = h_designDirection;
      medium_designDirection.Xi = Xi_designDirection;
      medium_nonDesignDirection.p = if flowDirection==Types.FlowDirection.Bidirectional then p_nonDesignDirection else Medium.p_default;
      medium_nonDesignDirection.h = if flowDirection==Types.FlowDirection.Bidirectional then h_nonDesignDirection else Medium.h_default;
      medium_nonDesignDirection.Xi = if flowDirection==Types.FlowDirection.Bidirectional then Xi_nonDesignDirection else zeros(Medium.nXi);
      
      assert(flowDirection==Types.FlowDirection.Bidirectional or m_flow>-0.001, "PartialTransportIsentropicAB: Mass flow direction was said to be positive only but is not.");
      
    end PartialTransportIsentropicAB;
    
    replaceable partial model PartialIdealJunction 
      "Partial infinitesimal junction model" 
      
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Fluid medium model" 
          annotation (choicesAllMatching=true);
      
      FluidPort_a port_1(redeclare package Medium = Medium) 
        annotation (extent=[-110,-10; -90,10]);
      FluidPort_a port_2(redeclare package Medium = Medium) 
        annotation (extent=[90,-10; 110,10]);
      FluidPort_a port_3(redeclare package Medium = Medium) 
        annotation (extent=[-10,90; 10,110]);
      
      // switch and declarations for manual tearing
      parameter Boolean useManualTearing = false 
        "Use manual tearing of ideal connection equations";
      
    end PartialIdealJunction;
    
    replaceable partial model PartialIdealJunctionAAB 
      "Partial infinitesimal junction model (two PortA's, one PortB, not supported for all interfaces)" 
      
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Fluid medium model" 
          annotation (choicesAllMatching=true);
      
      FluidPort_a port_1(redeclare package Medium = Medium) 
        annotation (extent=[-110,-10; -90,10]);
      FluidPort_b port_2(redeclare package Medium = Medium) 
        annotation (extent=[90,-10; 110,10]);
      FluidPort_a port_3(redeclare package Medium = Medium) 
        annotation (extent=[-10,90; 10,110]);
      
      // switch and declarations for manual tearing
      parameter Boolean useManualTearing = false 
        "Use manual tearing of ideal connection equations";
      
    end PartialIdealJunctionAAB;
    
    replaceable partial model PartialSource_A 
      "Partial source model with a Port_a" 
      
      // Medium model
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Medium model within the source" 
         annotation (choicesAllMatching=true);
      
      // Connector
      FluidPort_a port(redeclare package Medium = Medium) "Fluid connector" 
                        annotation (extent=[110,-10; 90,10]);
      
      // Medium instance  
      Medium.BaseProperties medium "Medium in the source";
      
      // Interface to generic implementation
    protected 
      Medium.MassFlowRate port_m_flow "Port mass flow rate";
      
      annotation (Icon(Text(extent=[-134,168; 134,106], string="%name")));
    end PartialSource_A;
    
    replaceable partial model PartialSource_B 
      "Partial source model with a Port_b" 
      
      // Medium model
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Medium model within the source" 
         annotation (choicesAllMatching=true);
      
      // Connector
      FluidPort_b port(redeclare package Medium = Medium) "Fluid connector" 
                        annotation (extent=[110,-10; 90,10]);
      
      // Medium instance  
      Medium.BaseProperties medium "Medium in the source";
      
      // Interface to generic implementation
    protected 
      Medium.MassFlowRate port_m_flow "Port mass flow rate";
      
      annotation (Icon(Text(extent=[-134,168; 134,106], string="%name")));
    end PartialSource_B;
    
    replaceable partial model PartialAsymmetricDistributedPipe 
      
      // Medium model
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Medium model in the component" 
         annotation (choicesAllMatching=true);
      
      //Fluid ports
      FluidPort_a port_a(redeclare package Medium = Medium) "Fluid inlet port" 
                       annotation (extent=[-110,-10; -90,10]);
      FluidPort_b port_b(redeclare package Medium = Medium) "Fluid outlet port"
                        annotation (extent=[90,-10; 110,10]);
      
      // Used for some approaches only
      SI.ThermalConductance G=0 
        "Conductive heat flow through component (used in some fluid interfaces only)";
      SI.DiffusionCoefficient H=0 
        "Mass diffusion through component (used in some fluid interfaces only)";
      
      annotation (Icon(Text(
            extent=[-148,-42; 148,-92],
            string="%name",
            style(gradient=2, fillColor=69))));
      
    end PartialAsymmetricDistributedPipe;
    
    replaceable partial model PartialSymmetricDistributedPipe 
      
      // Medium model
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Medium model in the component" 
         annotation (choicesAllMatching=true);
      
      //Fluid ports
      FluidPort_b port_a(redeclare package Medium = Medium) "Fluid inlet port" 
                       annotation (extent=[-110,-10; -90,10]);
      FluidPort_b port_b(redeclare package Medium = Medium) "Fluid outlet port"
                        annotation (extent=[90,-10; 110,10]);
      
      // Used for some approaches only
      SI.ThermalConductance G=0 
        "Conductive heat flow through component (used in some fluid interfaces only)";
      SI.DiffusionCoefficient H=0 
        "Mass diffusion through component (used in some fluid interfaces only)";
      
      annotation (Icon(Text(
            extent=[-148,-42; 148,-92],
            string="%name",
            style(gradient=2, fillColor=69))));
      
    end PartialSymmetricDistributedPipe;
    
    partial model PartialPortsAndMediumOnlyAB 
      
      // Medium model
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Medium model within the source" 
         annotation (choicesAllMatching=true);
      
      FluidPort_a port_a(redeclare package Medium = Medium) "Fluid inlet port" 
                         annotation (extent=[-112,-10; -92,10]);
      
      FluidPort_b port_b(redeclare package Medium = Medium) "Fluid outlet port"
                          annotation (extent=[90,-10; 110,10]);
    equation 
      
    end PartialPortsAndMediumOnlyAB;
    
    partial model PartialPortsAndMediumOnlyBA 
      
      // Medium model
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Medium model within the source" 
         annotation (choicesAllMatching=true);
      
      FluidPort_b port_a(redeclare package Medium = Medium) "Fluid inlet port" 
                         annotation (extent=[-112,-10; -92,10]);
      
      FluidPort_a port_b(redeclare package Medium = Medium) "Fluid outlet port"
                          annotation (extent=[90,-10; 110,10]);
    equation 
      
    end PartialPortsAndMediumOnlyBA;
  end PartialFluidInterface;
  
  partial package PartialFluidDiscretization 
    "Partial fluid discretization approach" 
    extends Icons.FluidDiscretization;
    
    constant Boolean isSymmetric 
      "Is this a symmetric discretization (with a half momentum balance at each end)?";
    
    replaceable partial model PartialDistributedFlow 
      
      // Nothing
      annotation (Icon(
             Rectangle(extent=[-100,40; 100,-40], style(
              color=69,
              gradient=2,
              fillColor=69))));
    end PartialDistributedFlow;
    
  end PartialFluidDiscretization;
  
  partial model PartialComponent 
    
    // Discretization
    replaceable package FluidDiscretization = PartialFluidDiscretization 
      "Fluid discretization used in this component" 
      annotation (choicesAllMatching = true);
    
    // Fluid interface
    replaceable package FluidInterface = PartialFluidInterface 
      "Fluid interface used in this component" 
      annotation (choicesAllMatching = true);
    
    // Medium model
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
      "Medium in the component"                                          annotation (
        choicesAllMatching =                                                                            true);
  end PartialComponent;
end Interfaces;
