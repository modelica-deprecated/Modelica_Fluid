package Common "Base classes common to more than one Component package" 
partial model PartialTwoPortTransport 
    "Partial element transporting fluid between two ports without storing mass or energy" 
    import SI = Modelica.SIunits;
    import Modelica.Constants;
  replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
                                                                       annotation (
      choicesAllMatching =                                                                            true);
    
 parameter Types.FlowDirectionWithGlobalDefault.Temp flowDirection=
                   Modelica_Fluid.Types.FlowDirection.Unidirectional 
      "Unidirectional (port_a -> port_b) or bidirectional flow component" 
     annotation(Dialog(tab="Advanced"));
    
  Interfaces.FluidPort_a port_a(redeclare package Medium = Medium,
                     m_flow(start=0,min=if allowFlowReversal then -Constants.inf else 0)) 
      "Fluid connector a (positive design flow direction is from port_a to port_b)"
    annotation (extent=[-110,-10; -90,10]);
  Interfaces.FluidPort_b port_b(redeclare package Medium = Medium,
                     m_flow(start=0,max=if allowFlowReversal then +Constants.inf else 0)) 
      "Fluid connector b (positive design flow direction is from port_a to port_b)"
    annotation (extent=[110,-10; 90,10]);
  Medium.BaseProperties medium_a "Medium properties in port_a";
  Medium.BaseProperties medium_b "Medium properties in port_b";
  Medium.MassFlowRate m_flow(start=0) 
      "Mass flow rate from port_a to port_b (m_flow > 0 is design flow direction)";
  SI.VolumeFlowRate V_flow_a = port_a.m_flow/medium_a.d 
      "Volume flow rate near port_a";
  SI.Pressure dp(start=0) "Pressure difference between port_a and port_b";
    
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
    
  /* Energy, mass and substance mass balance */
  port_a.H_flow + port_b.H_flow = 0;
  port_a.m_flow + port_b.m_flow = 0;
  port_a.mXi_flow + port_b.mXi_flow = zeros(Medium.nXi);
    
  // Design direction of mass flow rate
  m_flow = port_a.m_flow;
    
  // Pressure difference between ports
  dp = port_a.p - port_b.p;
    
end PartialTwoPortTransport;
  
  partial model PartialInitializationParameters 
    "Define parameter menu to initialize medium in component that has one medium model" 
    import Modelica_Fluid.Types;
    
    replaceable package Medium = PackageMedium extends 
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
    
    replaceable package Medium = PackageMedium extends 
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
    import Modelica.Constants;
    extends BaseClasses.Common.PartialInitializationParameters;
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching = true);
    
    parameter Types.FlowDirection.Temp flowDirection=
              Types.FlowDirection.Unidirectional 
      "Unidirectional (port_a -> port_b) or bidirectional flow component" 
       annotation(Dialog(tab="Advanced"));
    Interfaces.FluidPort_a port_a(redeclare package Medium = Medium,
                       m_flow(min=if allowFlowReversal then -Constants.inf else 0)) 
      "Fluid inlet port" annotation (extent=[-112,-10; -92,10]);
    Interfaces.FluidPort_b port_b(redeclare package Medium = Medium,
                       m_flow(max=if allowFlowReversal then +Constants.inf else 0)) 
      "Fluid outlet port" annotation (extent=[90,-10; 110,10]);
    Medium.BaseProperties medium(preferredMediumStates=true,
                 p(start=p_start), h(start=h_start),
                 T(start=T_start), Xi(start=X_start[1:Medium.nXi]));
    SI.Energy U "Internal energy of fluid";
    SI.Mass m "Mass of fluid";
    SI.Mass mXi[Medium.nXi] "Masses of independent components in the fluid";
    SI.Volume V_lumped "Volume";
    SI.HeatFlowRate Qs_flow "Heat flow across boundaries or energy source/sink";
    SI.Power Ws_flow "Work flow across boundaries or source term";
    annotation (
     Icon(
        Text(extent=[-144, 178; 146, 116], string="%name"),
        Text(
          extent=[-130, -108; 144, -150],
          style(color=0),
          string="V=%V")), Documentation(info="<html>
</html>"),
      Diagram);
    parameter Boolean allowFlowReversal=
       flowDirection == Types.FlowDirection.Bidirectional 
      "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
       annotation(Evaluate=true, Hide=true);
    equation 
    // boundary conditions
    port_a.p = medium.p;
    port_b.p = medium.p;
    port_a.H_flow = semiLinear(port_a.m_flow, port_a.h, medium.h);
    port_b.H_flow = semiLinear(port_b.m_flow, port_b.h, medium.h);
    port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, medium.Xi);
    port_b.mXi_flow = semiLinear(port_b.m_flow, port_b.Xi, medium.Xi);
    
    // Total quantities
    m    = V_lumped*medium.d;
    mXi = m*medium.Xi;
    U    = m*medium.u;
    
    // Mass and energy balance
    der(m)   = port_a.m_flow + port_b.m_flow;
    der(mXi) = port_a.mXi_flow + port_b.mXi_flow;
    der(U)   = port_a.H_flow + port_b.H_flow + Qs_flow +Ws_flow;
    
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
end Common;