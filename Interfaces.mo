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
              rgbfillColor={0,127,255})), Text(extent=[-150,110; 150,50], string=
                "%name")),
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
      
      import Modelica_Fluid.Types;
      
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
      Medium.BaseProperties medium(
        preferredMediumStates=true,
        p(start=p_start),
        h(start=h_start),
        T(start=T_start),
        Xi(start=X_start[1:Medium.nXi]));
      
    // Extensive properties
      SI.Energy U "Internal energy of fluid";
      SI.Mass m "Mass of fluid";
      SI.Mass mXi[Medium.nXi] "Masses of independent components in the fluid";
      SI.Volume V_lumped "Volume";
      
    //Initialization  
      parameter Types.Init.Temp initType=Types.Init.NoInit 
        "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p_start=Medium.p_default 
        "Start value of pressure" 
      annotation(Dialog(tab = "Initialization"));
      parameter Boolean use_T_start=true 
        "= true, use T_start, otherwise h_start" 
      annotation(Dialog(tab = "Initialization"), Evaluate=true);
      parameter Medium.Temperature T_start=if use_T_start then Medium.T_default else 
                Medium.temperature_phX(
            p_start,
            h_start,
            X_start) "Start value of temperature" 
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
      parameter Medium.SpecificEnthalpy h_start=if use_T_start then 
          Medium.specificEnthalpy_pTX(
            p_start,
            T_start,
            X_start) else Medium.h_default "Start value of specific enthalpy" 
      annotation(Dialog(tab = "Initialization", enable = not use_T_start));
      parameter Medium.MassFraction X_start[Medium.nX]=Medium.X_default 
        "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
      
      Medium.MassFlowRate m_flow_net "Net mass flow into the volume";
      Medium.MassFlowRate mXi_flow_net[Medium.nXi] 
        "Net substance mass flow into the volume";
      Medium.EnthalpyFlowRate H_flow_net "Net enthalpy flow into the volume";
      
    protected 
      SI.HeatFlowRate Qs_flow 
        "Heat flow across boundaries or energy source/sink";
      SI.Power Ws_flow "Work flow across boundaries or source term";
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
      
    equation 
    // Extensive quantities
      m = V_lumped*medium.d;
      mXi = m*medium.Xi;
      U = m*medium.u;
      
    // Mass and energy balance
      der(m) = m_flow_net;
      der(mXi) = mXi_flow_net;
      der(U) = H_flow_net + Qs_flow + Ws_flow;
      
    initial equation 
    // Initial conditions
      if initType == Types.Init.NoInit then
      // no initial equations
      elseif initType == Types.Init.InitialValues then
        if not Medium.singleState then
          medium.p = p_start;
        end if;
        if use_T_start then
          medium.T = T_start;
        else
          medium.h = h_start;
        end if;
        medium.Xi = X_start[1:Medium.nXi];
      elseif initType == Types.Init.SteadyState then
        if not Medium.singleState then
          der(medium.p) = 0;
        end if;
        der(medium.h) = 0;
        der(medium.Xi) = zeros(Medium.nXi);
      elseif initType == Types.Init.SteadyStateHydraulic then
        if not Medium.singleState then
          der(medium.p) = 0;
        end if;
        if use_T_start then
          medium.T = T_start;
        else
          medium.h = h_start;
        end if;
        medium.Xi = X_start[1:Medium.nXi];
      else
        assert(false, "Unsupported initialization option");
      end if;
    end PartialLumpedVolume;
    
    replaceable partial model PartialTwoPortTransport 
      "Partial isenthalpic element transporting fluid between two ports without storing mass or energy (two Port_b's)" 
      
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Medium in the component"                                          annotation (
          choicesAllMatching =                                                                            true);
      
      FluidPort_b port_a(redeclare package Medium = Medium) 
        "Fluid connector a (positive design flow direction is from port_a to port_b)"
        annotation (extent=[-110,-10; -90,10]);
      FluidPort_b port_b(redeclare package Medium = Medium) 
        "Fluid connector b (positive design flow direction is from port_a to port_b)"
        annotation (extent=[110,-10; 90,10]);
      
      Medium.MassFlowRate m_flow 
        "Mass flow rate from port_a to port_b (m_flow > 0 is design flow direction)";
      SI.Pressure dp "Pressure drop between port_a and port_b";
      
      // Used for some approaches only
      SI.ThermalConductance G=0 
        "Conductive heat flow through component (used in some fluid interfaces only)";
      SI.DiffusionCoefficient H=0 
        "Mass diffusion through component (used in some fluid interfaces only)";
      
      // Properties for fluid flowing in design direction, have to be provided by PartialTwoPortTransport
      // Used in pressure drop correlation
      Medium.AbsolutePressure p_designDirection 
        "Pressure for flow in design direction";
      Medium.SpecificEnthalpy h_designDirection 
        "Specific mixing enthalpy for flow in design direction";
      Medium.MassFraction Xi_designDirection[Medium.nXi] 
        "Mass fractions for flow in design direction";
      
      // Properties for fluid flowing in non-design direction, have to be provided by PartialTwoPortTransport
      // Used in pressure drop correlation
      Medium.AbsolutePressure p_nonDesignDirection 
        "Pressure for flow in non-design direction";
      Medium.SpecificEnthalpy h_nonDesignDirection 
        "Specific mixing enthalpy for flow in non-design direction";
      Medium.MassFraction Xi_nonDesignDirection[Medium.nXi] 
        "Mass fractions for flow in non-design direction";
      
      annotation (Icon(Text(
            extent=[-144,119; 144,70],
            string="%name",
            style(gradient=2, fillColor=69))));
    equation 
      
    end PartialTwoPortTransport;
    
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
    end PartialIdealJunction;
    
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

    replaceable partial model PartialTwoPortTransportAA 
      "Partial isenthalpic element transporting fluid between two ports without storing mass or energy (two Port_a's, not supported for all interfaces)" 
      
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Medium in the component"                                          annotation (
          choicesAllMatching =                                                                            true);
      
      FluidPort_a port_a(redeclare package Medium = Medium) 
        "Fluid connector a (positive design flow direction is from port_a to port_b)"
        annotation (extent=[-110,-10; -90,10]);
      FluidPort_a port_b(redeclare package Medium = Medium) 
        "Fluid connector b (positive design flow direction is from port_a to port_b)"
        annotation (extent=[110,-10; 90,10]);
      
      Medium.MassFlowRate m_flow 
        "Mass flow rate from port_a to port_b (m_flow > 0 is design flow direction)";
      SI.Pressure dp "Pressure drop between port_a and port_b";
      
      // Used for some approaches only
      SI.ThermalConductance G=0 
        "Conductive heat flow through component (used in some fluid interfaces only)";
      SI.DiffusionCoefficient H=0 
        "Mass diffusion through component (used in some fluid interfaces only)";
      
      // Properties for fluid flowing in design direction, have to be provided by PartialTwoPortTransport
      // Used in pressure drop correlation
      Medium.AbsolutePressure p_designDirection 
        "Pressure for flow in design direction";
      Medium.SpecificEnthalpy h_designDirection 
        "Specific mixing enthalpy for flow in design direction";
      Medium.MassFraction Xi_designDirection[Medium.nXi] 
        "Mass fractions for flow in design direction";
      
      // Properties for fluid flowing in non-design direction, have to be provided by PartialTwoPortTransport
      // Used in pressure drop correlation
      Medium.AbsolutePressure p_nonDesignDirection 
        "Pressure for flow in non-design direction";
      Medium.SpecificEnthalpy h_nonDesignDirection 
        "Specific mixing enthalpy for flow in non-design direction";
      Medium.MassFraction Xi_nonDesignDirection[Medium.nXi] 
        "Mass fractions for flow in non-design direction";
      
      annotation (Icon(Text(
            extent=[-144,119; 144,70],
            string="%name",
            style(gradient=2, fillColor=69))));
    equation 
      
    end PartialTwoPortTransportAA;
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
