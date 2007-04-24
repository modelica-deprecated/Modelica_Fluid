package Interfaces 
  "Interfaces for steady state and unsteady, mixed-phase, multi-substance, incompressible and compressible flow" 
  
  annotation (Documentation(info="<html>
 
</html>", revisions="<html>
<ul>
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
  
partial model PartialTwoPortTransport 
    "Partial element transporting fluid between two ports without storing mass or energy" 
      import Modelica.Constants;
  replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
                                                                       annotation (
      choicesAllMatching =                                                                            true);
    
  //Initialization
  parameter Medium.AbsolutePressure p_a_start "Guess value of pressure" 
    annotation(Dialog(tab = "Guess Value Initialization"));
  parameter Medium.AbsolutePressure p_b_start "Guess value of pressure" 
    annotation(Dialog(tab = "Guess Value Initialization"));
  parameter Boolean use_T_start = true "= true, use T_start, otherwise h_start"
    annotation(Dialog(tab = "Guess Value Initialization"), Evaluate=true);
  parameter Medium.Temperature T_start=
    if use_T_start then Medium.T_default else Medium.temperature_phX(p_a_start,h_start,X_start) 
      "Guess value of temperature" 
    annotation(Dialog(tab = "Guess Value Initialization", enable = use_T_start));
  parameter Medium.SpecificEnthalpy h_start=
    if use_T_start then Medium.specificEnthalpy_pTX(p_a_start, T_start, X_start) else Medium.h_default 
      "Guess value of specific enthalpy" 
    annotation(Dialog(tab = "Guess Value Initialization", enable = not use_T_start));
  parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
      "Guess value of mass fractions m_i/m" 
    annotation (Dialog(tab="Guess Value Initialization", enable=Medium.nXi > 0));
 parameter Types.FlowDirection.Temp flowDirection=
                   Modelica_Fluid.Types.FlowDirection.Bidirectional 
      "Unidirectional (port_a -> port_b) or bidirectional flow component" 
     annotation(Dialog(tab="Advanced"));
    
  Modelica.Fluid.Interfaces.FluidPort_a port_a(
                                redeclare package Medium = Medium,
                     m_flow(start=0,min=if allowFlowReversal then -Constants.inf else 0)) 
      "Fluid connector a (positive design flow direction is from port_a to port_b)"
    annotation (extent=[-110,-10; -90,10]);
  Modelica.Fluid.Interfaces.FluidPort_b port_b(
                                redeclare package Medium = Medium,
                     m_flow(start=0,max=if allowFlowReversal then +Constants.inf else 0)) 
      "Fluid connector b (positive design flow direction is from port_a to port_b)"
    annotation (extent=[110,-10; 90,10]);
  Medium.BaseProperties medium_a(p(start=p_a_start), h(start=h_start), X(start=X_start)) 
      "Medium properties in port_a";
  Medium.BaseProperties medium_b(p(start=p_b_start), h(start=h_start), X(start=X_start)) 
      "Medium properties in port_b";
  Medium.MassFlowRate m_flow 
      "Mass flow rate from port_a to port_b (m_flow > 0 is design flow direction)";
  SI.VolumeFlowRate V_flow_a = port_a.m_flow/medium_a.d 
      "Volume flow rate near port_a";
  SI.Pressure dp(start=p_a_start-p_b_start) 
      "Pressure difference between port_a and port_b";
    
  annotation (
    Coordsys(grid=[1, 1], component=[20, 20]),
    Diagram,
    Documentation(info="<html>
<p>
This component transports fluid between its two ports, without
storing mass or energy. Reversal and zero mass flow rate is taken
care of, for details see definition of built-in operator semiLinear().
<p>
When using this partial component, an equation for the momentum
balance has to be added by specifying a relationship
between the pressure drop <tt>dp</tt> and the mass flow rate <tt>m_flow</tt>.
</p>
</html>"),
    Icon);
  protected 
    parameter Boolean allowFlowReversal=
     flowDirection == Modelica_Fluid.Types.FlowDirection.Bidirectional 
      "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
     annotation(Evaluate=true, Hide=true);
equation 
  // Properties in the ports
  port_a.p   = medium_a.p;
  port_a.h   = medium_a.h;
  port_a.Xi = medium_a.Xi;
  port_b.p   = medium_b.p;
  port_b.h   = medium_b.h;
  port_b.Xi = medium_b.Xi;
    
  /* Handle reverse and zero flow */
  port_a.H_flow   = semiLinear(port_a.m_flow, port_a.h,  port_b.h);
  port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, port_b.Xi);
  port_a.mC_flow = semiLinear(port_a.m_flow, port_a.C, port_b.C);
    
  /* Energy, mass and substance mass balance */
  port_a.H_flow + port_b.H_flow = 0;
  port_a.m_flow + port_b.m_flow = 0;
  port_a.mXi_flow + port_b.mXi_flow = zeros(Medium.nXi);
  port_a.mC_flow + port_b.mC_flow = zeros(Medium.nC);
    
  // Design direction of mass flow rate
  m_flow = port_a.m_flow;
    
  // Pressure difference between ports
  dp = port_a.p - port_b.p;
    
end PartialTwoPortTransport;
  
  partial model PartialInitializationParameters 
    "Define parameter menu to initialize medium in component that has one medium model" 
    import Modelica_Fluid.Types;
    
    replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching = true);
    parameter Types.Init.Temp initType=
              Types.Init.NoInit "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_start = Medium.p_default 
      "Start value of pressure" 
      annotation(Dialog(tab = "Initialization"));
    parameter Boolean use_T_start = true 
      "= true, use T_start, otherwise h_start" 
      annotation(Dialog(tab = "Initialization"), Evaluate=true);
    parameter Medium.Temperature T_start=
      if use_T_start then Medium.T_default else Medium.temperature_phX(p_start,h_start,X_start) 
      "Start value of temperature" 
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.specificEnthalpy_pTX(p_start, T_start, X_start) else Medium.h_default 
      "Start value of specific enthalpy" 
      annotation(Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
  end PartialInitializationParameters;
  
  partial model PartialGuessValueParameters 
    "Define parameter menu to initialize guess values of medium in component that has one medium model" 
    import Modelica_Fluid.Types;
    
    replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching = true);
    parameter Medium.AbsolutePressure p_start = Medium.p_default 
      "Guess value of pressure" 
      annotation(Dialog(tab = "Guess Value Initialization"));
    parameter Boolean use_T_start = true 
      "= true, use T_start, otherwise h_start" 
      annotation(Dialog(tab = "Guess Value Initialization"), Evaluate=true);
    parameter Medium.Temperature T_start=
      if use_T_start then Medium.T_default else Medium.temperature_phX(p_start,h_start,X_start) 
      "Guess value of temperature" 
      annotation(Dialog(tab = "Guess Value Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.specificEnthalpy_pTX(p_start, T_start, X_start) else Medium.h_default 
      "Guess value of specific enthalpy" 
      annotation(Dialog(tab = "Guess Value Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
      "Guess value of mass fractions m_i/m" 
      annotation (Dialog(tab="Guess Value Initialization", enable=Medium.nXi > 0));
  end PartialGuessValueParameters;
  
    partial model PartialLumpedVolume 
    "Mixing volume with inlet and outlet ports (flow reversal is allowed)" 
      extends Modelica_Fluid.Interfaces.PartialInitializationParameters;
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
      "Medium in the component" 
          annotation (choicesAllMatching = true);
    
      parameter Types.FlowDirection.Temp flowDirection=Types.FlowDirection.
          Bidirectional 
      "Unidirectional (port_a -> port_b) or bidirectional flow component" 
       annotation(Dialog(tab="Advanced"));
       parameter Boolean allowFlowReversal=flowDirection == Types.FlowDirection.
          Bidirectional 
      "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
       annotation(Evaluate=true, Hide=true);
    
      Modelica.Fluid.Interfaces.FluidPort_a port_a(
                                    redeclare package Medium = Medium, m_flow(min=
              if allowFlowReversal then -Modelica.Constants.inf else 0)) 
      "Fluid inlet port" annotation (extent=[-112,-10; -92,10]);
      Modelica.Fluid.Interfaces.FluidPort_b port_b(
                                    redeclare package Medium = Medium, m_flow(max=
              if allowFlowReversal then +Modelica.Constants.inf else 0)) 
      "Fluid outlet port" annotation (extent=[90,-10; 110,10]);
      Medium.BaseProperties medium(
        preferredMediumStates=true,
        p(start=p_start),
        h(start=h_start),
        T(start=T_start),
        Xi(start=X_start[1:Medium.nXi]));
      SI.Energy U "Internal energy of fluid";
      SI.Mass m "Mass of fluid";
      SI.Mass mXi[Medium.nXi] "Masses of independent components in the fluid";
      SI.Volume V_lumped "Volume";
    
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
    // boundary conditions
      port_a.p = medium.p;
      port_b.p = medium.p;
      port_a.H_flow = semiLinear(
        port_a.m_flow,
        port_a.h,
        medium.h);
      port_b.H_flow = semiLinear(
        port_b.m_flow,
        port_b.h,
        medium.h);
      port_a.mXi_flow = semiLinear(
        port_a.m_flow,
        port_a.Xi,
        medium.Xi);
      port_a.mC_flow = semiLinear(
        port_a.m_flow,
        port_a.C,
        port_b.C);
      port_b.mXi_flow = semiLinear(
        port_b.m_flow,
        port_b.Xi,
        medium.Xi);
    
    // Total quantities
      m = V_lumped*medium.d;
      mXi = m*medium.Xi;
      U = m*medium.u;
    
    // Mass and energy balance
      der(m) = port_a.m_flow + port_b.m_flow;
      der(mXi) = port_a.mXi_flow + port_b.mXi_flow;
      der(U) = port_a.H_flow + port_b.H_flow + Qs_flow + Ws_flow;
      zeros(Medium.nC) = port_a.mC_flow+port_b.mC_flow;
    
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
end Interfaces;
