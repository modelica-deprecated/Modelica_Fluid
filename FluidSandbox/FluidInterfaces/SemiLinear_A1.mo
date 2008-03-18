within FluidSandbox.FluidInterfaces;
package SemiLinear_A1 
  "Implementation A1 using SemiLinear (similar to the Modelica_Fluid library)" 
  extends Interfaces.PartialFluidInterface(usesNewConnectionSemantics=false);
  redeclare replaceable connector extends FluidPort 
    "Interface for quasi one-dimensional fluid flow in a piping network (incompressible or compressible, one or more phases, one or more substances)" 
    
    Medium.AbsolutePressure p "Pressure in the connection point";
    flow Medium.MassFlowRate m_flow 
      "Mass flow rate from the connection point into the component";
    
    Medium.SpecificEnthalpy h 
      "Specific mixing enthalpy in the connection point";
    flow Medium.EnthalpyFlowRate H_flow 
      "Enthalpy flow rate into the component (if m_flow > 0, H_flow = m_flow*h)";
    
    Medium.MassFraction Xi[Medium.nXi] 
      "Independent mixture mass fractions m_i/m in the connection point";
    flow Medium.MassFlowRate mXi_flow[Medium.nXi] 
      "Mass flow rates of the independent substances from the connection point into the component (if m_flow > 0, mXi_flow = m_flow*Xi)";
    
  end FluidPort;
  
  redeclare replaceable connector extends FluidPort_a 
    "Generic fluid connector at design inlet" 
    
    // Inheritance of FluidPort is sufficient
  end FluidPort_a;
  
  redeclare replaceable connector extends FluidPort_b 
    "Generic fluid connector at design outlet" 
    
    // Inheritance of FluidPort is sufficient
  end FluidPort_b;
  
  redeclare replaceable model extends ConnectionSemantics 
    "No unsupported connection semantics for this fluid interface" 
    
  equation 
    // No unsupported connection semantics for this approach
    
    annotation (
      defaultComponentName="semantics",
      Coordsys(extent=[-100,-20; 100,20], scale=0.05),
      Icon(Text(extent=[-300,0; 300,-64], string="unused")),
      Diagram);
  end ConnectionSemantics;
  
  redeclare replaceable partial model extends PartialLumpedVolume 
    "Mixing volume with inlet and outlet ports (flow reversal is allowed)" 
    
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
  // Pressure
    port_a.p = medium.p*ones(n_a);
    port_b.p = medium.p*ones(n_b);
    
  // Enthalpy flow
    port_a.H_flow = semiLinear(
            port_a.m_flow,
            port_a.h,
            h_a_outflow*ones(n_a));
    port_b.H_flow = semiLinear(
            port_b.m_flow,
            port_b.h,
            h_b_outflow*ones(n_b));
    
  // Substance mass flows
    port_a.mXi_flow = {semiLinear(
            port_a[i].m_flow,
            port_a[i].Xi,
            Xi_a_outflow) for i in 1:n_a};
    port_b.mXi_flow = {semiLinear(
            port_b[i].m_flow,
            port_b[i].Xi,
            Xi_b_outflow) for i in 1:n_b};
      /*
  for i in 1:n_a loop
  port_a[i].mXi_flow = semiLinear(
    port_a[i].m_flow,
    port_a[i].Xi,
    medium.Xi);
  end for;
  for i in 1:n_a loop
  port_b[i].mXi_flow = semiLinear(
    port_b[i].m_flow,
    port_b[i].Xi,
    medium.Xi);
  end for;
*/
    
  // Net flow rates
    m_flow_net = sum(port_a.m_flow) + sum(port_b.m_flow);
    mXi_flow_net = {sum(port_a[:].mXi_flow[i]) + sum(port_b[:].mXi_flow[i]) 
      for i in 1:Medium.nXi};
    H_flow_net = sum(port_a.H_flow) + sum(port_b.H_flow);
      /*
  Each substance separately
  for i in 1:Medium.nXi loop
    mXi_flow_net[i] = sum(port_a[:].mXi_flow[i]) + sum(port_b[:].mXi_flow[i]);
  end for;
*/
    
  end PartialLumpedVolume;
  
  
  redeclare replaceable partial model extends PartialTransportIsenthalpic 
    "Partial isenthalpic element transporting fluid between two ports without storing mass or energy (two Port_b's)" 
    
  equation 
  /* Handle reverse and zero flow */
    port_a.H_flow = semiLinear(
            port_a.m_flow,
            port_a.h,
            port_b.h);
    port_a.mXi_flow = semiLinear(
            port_a.m_flow,
            port_a.Xi,
            port_b.Xi);
    
  /* Mass, energy, substance mass balance */
    port_a.m_flow + port_b.m_flow = 0;
    port_a.H_flow + port_b.H_flow = 0;
    port_a.mXi_flow + port_b.mXi_flow = zeros(Medium.nXi);
    
  // Design direction of mass flow rate
    m_flow = port_a.m_flow;
    
  // Pressure difference between ports
    dp = port_a.p - port_b.p;
    
    // This approach does not provide both potential upstream properties 
    // (only upstream properties for current mass flow direction and downstream mixing properties)
    p_designDirection = port_a.p 
      "Upstream pressure if flow was in design direction";
    h_designDirection = port_a.h 
      "Approximation of upstream specific enthalpy if flow was in design direction";
    Xi_designDirection = port_a.Xi 
      "Approxiamtion of upstream mass fractions if flow was in design direction";
    p_nonDesignDirection = port_b.p 
      "Upstream pressure if flow was in non-design direction";
    h_nonDesignDirection = port_b.h 
      "Approximation of upstream specific enthalpy if flow was in non-design direction";
    Xi_nonDesignDirection = port_b.Xi 
      "Approximation of upstream mass fractions if flow was in non-design direction";
    
    // sensors
    calc_T_a = if provide_T_a then calc_T_a_medium.T else 0;
    calc_T_b = if provide_T_b then calc_T_b_medium.T else 0;
    calc_p_a = if provide_p_a then port_a.p else 0;
    calc_p_b = if provide_p_b then port_b.p else 0;
    calc_m_flow_ab = if provide_m_flow_ab then m_flow else 0;
    
    calc_T_a_medium.p = if provide_T_a then port_a.p else Medium.p_default;
    calc_T_a_medium.h = if provide_T_a then (if port_a.m_flow > 0 then port_a.h else port_b.h) else Medium.h_default;
    calc_T_a_medium.Xi = if provide_T_a then (if port_a.m_flow > 0 then port_a.Xi else port_b.Xi) else zeros(Medium.nXi);
    calc_T_b_medium.p = if provide_T_b then port_b.p else Medium.p_default;
    calc_T_b_medium.h = if provide_T_b then (if port_b.m_flow > 0 then port_b.h else port_a.h) else Medium.h_default;
    calc_T_b_medium.Xi = if provide_T_b then (if port_b.m_flow > 0 then port_b.Xi else port_a.Xi) else zeros(Medium.nXi);
    /*
  // Simpler approach for temperature like that below results in numeric Jacobian
  calc_T_a = if provide_T_a then 
    Medium.temperature(Medium.setState_phX(
      port_a.p, if port_a.m_flow > 0 then port_a.h else port_b.h,
      if port_a.m_flow > 0 then port_a.Xi else port_b.Xi)) else 
         0 "Results in numeric Jacobian";
  calc_T_b = if provide_T_b then 
    Medium.temperature(Medium.setState_phX(
      port_b.p, if port_b.m_flow > 0 then port_b.h else port_a.h,
      if port_b.m_flow > 0 then port_b.Xi else port_a.Xi)) else 
         0 "Results in numeric Jacobian";
  */
  end PartialTransportIsenthalpic;
  
  redeclare replaceable partial model extends PartialTransportIsenthalpicAA 
    "Partial isenthalpic element transporting fluid between two ports without storing mass or energy (two Port_a's, allowed in this approach)" 
    
  equation 
  /* Handle reverse and zero flow */
    port_a.H_flow = semiLinear(
            port_a.m_flow,
            port_a.h,
            port_b.h);
    port_a.mXi_flow = semiLinear(
            port_a.m_flow,
            port_a.Xi,
            port_b.Xi);
    
  /* Mass, energy, substance mass balance */
    port_a.m_flow + port_b.m_flow = 0;
    port_a.H_flow + port_b.H_flow = 0;
    port_a.mXi_flow + port_b.mXi_flow = zeros(Medium.nXi);
    
  // Design direction of mass flow rate
    m_flow = port_a.m_flow;
    
  // Pressure difference between ports
    dp = port_a.p - port_b.p;
    
    // This approach does not provide both potential upstream properties 
    // (only upstream properties for current mass flow direction and downstream mixing properties)
    p_designDirection = port_a.p 
      "Upstream pressure if flow was in design direction";
    h_designDirection = port_a.h 
      "Approximation of upstream specific enthalpy if flow was in design direction";
    Xi_designDirection = port_a.Xi 
      "Approxiamtion of upstream mass fractions if flow was in design direction";
    p_nonDesignDirection = port_b.p 
      "Upstream pressure if flow was in non-design direction";
    h_nonDesignDirection = port_b.h 
      "Approximation of upstream specific enthalpy if flow was in non-design direction";
    Xi_nonDesignDirection = port_b.Xi 
      "Approximation of upstream mass fractions if flow was in non-design direction";
    
    // sensors
    calc_T_a = if provide_T_a then calc_T_a_medium.T else 0;
    calc_T_b = if provide_T_b then calc_T_b_medium.T else 0;
    calc_p_a = if provide_p_a then port_a.p else 0;
    calc_p_b = if provide_p_b then port_b.p else 0;
    calc_m_flow_ab = if provide_m_flow_ab then m_flow else 0;
    
    calc_T_a_medium.p = if provide_T_a then port_a.p else Medium.p_default;
    calc_T_a_medium.h = if provide_T_a then (if port_a.m_flow > 0 then port_a.h else port_b.h) else Medium.h_default;
    calc_T_a_medium.Xi = if provide_T_a then (if port_a.m_flow > 0 then port_a.Xi else port_b.Xi) else zeros(Medium.nXi);
    calc_T_b_medium.p = if provide_T_b then port_b.p else Medium.p_default;
    calc_T_b_medium.h = if provide_T_b then (if port_b.m_flow > 0 then port_b.h else port_a.h) else Medium.h_default;
    calc_T_b_medium.Xi = if provide_T_b then (if port_b.m_flow > 0 then port_b.Xi else port_a.Xi) else zeros(Medium.nXi);
    /*
  // Simpler approach for temperature like that below results in numeric Jacobian
  calc_T_a = if provide_T_a then 
    Medium.temperature(Medium.setState_phX(
      port_a.p, if port_a.m_flow > 0 then port_a.h else port_b.h,
      if port_a.m_flow > 0 then port_a.Xi else port_b.Xi)) else 
         0 "Results in numeric Jacobian";
  calc_T_b = if provide_T_b then 
    Medium.temperature(Medium.setState_phX(
      port_b.p, if port_b.m_flow > 0 then port_b.h else port_a.h,
      if port_b.m_flow > 0 then port_b.Xi else port_a.Xi)) else 
         0 "Results in numeric Jacobian";
  */
  end PartialTransportIsenthalpicAA;
  
  redeclare replaceable partial model extends PartialTransportIsenthalpicAB 
    "Partial isenthalpic element transporting fluid between two ports without storing mass or energy (a Port_a and Port_b each, allowed in this approach)" 
    
  equation 
  /* Handle reverse and zero flow */
    port_a.H_flow = semiLinear(
            port_a.m_flow,
            port_a.h,
            port_b.h);
    port_a.mXi_flow = semiLinear(
            port_a.m_flow,
            port_a.Xi,
            port_b.Xi);
    
  /* Mass, energy, substance mass balance */
    port_a.m_flow + port_b.m_flow = 0;
    port_a.H_flow + port_b.H_flow = 0;
    port_a.mXi_flow + port_b.mXi_flow = zeros(Medium.nXi);
    
  // Design direction of mass flow rate
    m_flow = port_a.m_flow;
    
  // Pressure difference between ports
    dp = port_a.p - port_b.p;
    
    // This approach does not provide both potential upstream properties 
    // (only upstream properties for current mass flow direction and downstream mixing properties)
    p_designDirection = port_a.p 
      "Upstream pressure if flow was in design direction";
    h_designDirection = port_a.h 
      "Approximation of upstream specific enthalpy if flow was in design direction";
    Xi_designDirection = port_a.Xi 
      "Approxiamtion of upstream mass fractions if flow was in design direction";
    p_nonDesignDirection = port_b.p 
      "Upstream pressure if flow was in non-design direction";
    h_nonDesignDirection = port_b.h 
      "Approximation of upstream specific enthalpy if flow was in non-design direction";
    Xi_nonDesignDirection = port_b.Xi 
      "Approximation of upstream mass fractions if flow was in non-design direction";
    
    // sensors
    calc_T_a = if provide_T_a then calc_T_a_medium.T else 0;
    calc_T_b = if provide_T_b then calc_T_b_medium.T else 0;
    calc_p_a = if provide_p_a then port_a.p else 0;
    calc_p_b = if provide_p_b then port_b.p else 0;
    calc_m_flow_ab = if provide_m_flow_ab then m_flow else 0;
    
    calc_T_a_medium.p = if provide_T_a then port_a.p else Medium.p_default;
    calc_T_a_medium.h = if provide_T_a then (if port_a.m_flow > 0 then port_a.h else port_b.h) else Medium.h_default;
    calc_T_a_medium.Xi = if provide_T_a then (if port_a.m_flow > 0 then port_a.Xi else port_b.Xi) else zeros(Medium.nXi);
    calc_T_b_medium.p = if provide_T_b then port_b.p else Medium.p_default;
    calc_T_b_medium.h = if provide_T_b then (if port_b.m_flow > 0 then port_b.h else port_a.h) else Medium.h_default;
    calc_T_b_medium.Xi = if provide_T_b then (if port_b.m_flow > 0 then port_b.Xi else port_a.Xi) else zeros(Medium.nXi);
    /*
  // Simpler approach for temperature like that below results in numeric Jacobian
  calc_T_a = if provide_T_a then 
    Medium.temperature(Medium.setState_phX(
      port_a.p, if port_a.m_flow > 0 then port_a.h else port_b.h,
      if port_a.m_flow > 0 then port_a.Xi else port_b.Xi)) else 
         0 "Results in numeric Jacobian";
  calc_T_b = if provide_T_b then 
    Medium.temperature(Medium.setState_phX(
      port_b.p, if port_b.m_flow > 0 then port_b.h else port_a.h,
      if port_b.m_flow > 0 then port_b.Xi else port_a.Xi)) else 
         0 "Results in numeric Jacobian";
  */
  end PartialTransportIsenthalpicAB;
  
  redeclare replaceable partial model extends PartialTransportIsentropic 
    "Partial isentropic element transporting fluid between two ports without storing mass or energy (two Port_b's)" 
    
    import Modelica_Fluid.Types;
    
    Medium.SpecificEnthalpy h_a_outflow = if flowDirection==Types.FlowDirection.Bidirectional then port_b.h - eta_ise*(port_b.h - Medium.isentropicEnthalpy(port_a.p, medium_nonDesignDirection.state)) else port_b.h;
    Medium.SpecificEnthalpy h_b_outflow = port_a.h - eta_ise*(port_a.h - Medium.isentropicEnthalpy(port_b.p, medium_designDirection.state));
    
  equation 
  /* Handle reverse and zero flow */
    port_a.H_flow = semiLinear(
            port_a.m_flow,
            port_a.h,
            h_a_outflow);
    port_b.H_flow = semiLinear(port_b.m_flow,
            port_b.h,
            h_b_outflow);
    port_a.mXi_flow = semiLinear(
            port_a.m_flow,
            port_a.Xi,
            port_b.Xi);
    
  /* Mass, energy, substance mass balance */
    port_a.m_flow + port_b.m_flow = 0;
    port_a.H_flow + port_b.H_flow + P_mechanical = 0;
    port_a.mXi_flow + port_b.mXi_flow = zeros(Medium.nXi);
    
  // Design direction of mass flow rate
    m_flow = port_a.m_flow;
    
  // Pressure difference between ports
    dp = port_a.p - port_b.p;
    
    // This approach does not provide both potential upstream properties 
    // (only upstream properties for current mass flow direction and downstream mixing properties)
    p_designDirection = port_a.p 
      "Upstream pressure if flow was in design direction";
    h_designDirection = port_a.h 
      "Approximation of upstream specific enthalpy if flow was in design direction";
    Xi_designDirection = port_a.Xi 
      "Approxiamtion of upstream mass fractions if flow was in design direction";
    p_nonDesignDirection = port_b.p 
      "Upstream pressure if flow was in non-design direction";
    h_nonDesignDirection = port_b.h 
      "Approximation of upstream specific enthalpy if flow was in non-design direction";
    Xi_nonDesignDirection = port_b.Xi 
      "Approximation of upstream mass fractions if flow was in non-design direction";
    
    // sensors
    calc_T_a = if provide_T_a then calc_T_a_medium.T else 0;
    calc_T_b = if provide_T_b then calc_T_b_medium.T else 0;
    calc_p_a = if provide_p_a then port_a.p else 0;
    calc_p_b = if provide_p_b then port_b.p else 0;
    calc_m_flow_ab = if provide_m_flow_ab then m_flow else 0;
    
    calc_T_a_medium.p = if provide_T_a then port_a.p else Medium.p_default;
    calc_T_a_medium.h = if provide_T_a then (if port_a.m_flow > 0 then port_a.h else h_a_outflow) else Medium.h_default;
    calc_T_a_medium.Xi = if provide_T_a then (if port_a.m_flow > 0 then port_a.Xi else port_b.Xi) else zeros(Medium.nXi);
    calc_T_b_medium.p = if provide_T_b then port_b.p else Medium.p_default;
    calc_T_b_medium.h = if provide_T_b then (if port_b.m_flow > 0 then port_b.h else h_b_outflow) else Medium.h_default;
    calc_T_b_medium.Xi = if provide_T_b then (if port_b.m_flow > 0 then port_b.Xi else port_a.Xi) else zeros(Medium.nXi);
    /*
  // Simpler approach for temperature like that below results in numeric Jacobian
  calc_T_a = if provide_T_a then 
    Medium.temperature(Medium.setState_phX(
      port_a.p, if port_a.m_flow > 0 then port_a.h else port_b.h,
      if port_a.m_flow > 0 then port_a.Xi else port_b.Xi)) else 
         0 "Results in numeric Jacobian";
  calc_T_b = if provide_T_b then 
    Medium.temperature(Medium.setState_phX(
      port_b.p, if port_b.m_flow > 0 then port_b.h else port_a.h,
      if port_b.m_flow > 0 then port_b.Xi else port_a.Xi)) else 
         0 "Results in numeric Jacobian";
  */
    
  end PartialTransportIsentropic;
  
  redeclare replaceable partial model extends PartialTransportIsentropicAA 
    "Partial isentropic element transporting fluid between two ports without storing mass or energy (two Port_a's, allowed in this approach)" 
    
    import Modelica_Fluid.Types;
    
    Medium.SpecificEnthalpy h_a_outflow = if flowDirection==Types.FlowDirection.Bidirectional then port_b.h - eta_ise*(port_b.h - Medium.isentropicEnthalpy(port_a.p, medium_nonDesignDirection.state)) else port_b.h;
    Medium.SpecificEnthalpy h_b_outflow = port_a.h - eta_ise*(port_a.h - Medium.isentropicEnthalpy(port_b.p, medium_designDirection.state));
    
  equation 
  /* Handle reverse and zero flow */
    port_a.H_flow = semiLinear(
            port_a.m_flow,
            port_a.h,
            h_a_outflow);
    port_b.H_flow = semiLinear(port_b.m_flow,
            port_b.h,
            h_b_outflow);
    port_a.mXi_flow = semiLinear(
            port_a.m_flow,
            port_a.Xi,
            port_b.Xi);
    
  /* Mass, energy, substance mass balance */
    port_a.m_flow + port_b.m_flow = 0;
    port_a.H_flow + port_b.H_flow + P_mechanical = 0;
    port_a.mXi_flow + port_b.mXi_flow = zeros(Medium.nXi);
    
  // Design direction of mass flow rate
    m_flow = port_a.m_flow;
    
  // Pressure difference between ports
    dp = port_a.p - port_b.p;
    
    // This approach does not provide both potential upstream properties 
    // (only upstream properties for current mass flow direction and downstream mixing properties)
    p_designDirection = port_a.p 
      "Upstream pressure if flow was in design direction";
    h_designDirection = port_a.h 
      "Approximation of upstream specific enthalpy if flow was in design direction";
    Xi_designDirection = port_a.Xi 
      "Approxiamtion of upstream mass fractions if flow was in design direction";
    p_nonDesignDirection = port_b.p 
      "Upstream pressure if flow was in non-design direction";
    h_nonDesignDirection = port_b.h 
      "Approximation of upstream specific enthalpy if flow was in non-design direction";
    Xi_nonDesignDirection = port_b.Xi 
      "Approximation of upstream mass fractions if flow was in non-design direction";
    
    // sensors
    calc_T_a = if provide_T_a then calc_T_a_medium.T else 0;
    calc_T_b = if provide_T_b then calc_T_b_medium.T else 0;
    calc_p_a = if provide_p_a then port_a.p else 0;
    calc_p_b = if provide_p_b then port_b.p else 0;
    calc_m_flow_ab = if provide_m_flow_ab then m_flow else 0;
    
    calc_T_a_medium.p = if provide_T_a then port_a.p else Medium.p_default;
    calc_T_a_medium.h = if provide_T_a then (if port_a.m_flow > 0 then port_a.h else h_a_outflow) else Medium.h_default;
    calc_T_a_medium.Xi = if provide_T_a then (if port_a.m_flow > 0 then port_a.Xi else port_b.Xi) else zeros(Medium.nXi);
    calc_T_b_medium.p = if provide_T_b then port_b.p else Medium.p_default;
    calc_T_b_medium.h = if provide_T_b then (if port_b.m_flow > 0 then port_b.h else h_b_outflow) else Medium.h_default;
    calc_T_b_medium.Xi = if provide_T_b then (if port_b.m_flow > 0 then port_b.Xi else port_a.Xi) else zeros(Medium.nXi);
    /*
  // Simpler approach for temperature like that below results in numeric Jacobian
  calc_T_a = if provide_T_a then 
    Medium.temperature(Medium.setState_phX(
      port_a.p, if port_a.m_flow > 0 then port_a.h else port_b.h,
      if port_a.m_flow > 0 then port_a.Xi else port_b.Xi)) else 
         0 "Results in numeric Jacobian";
  calc_T_b = if provide_T_b then 
    Medium.temperature(Medium.setState_phX(
      port_b.p, if port_b.m_flow > 0 then port_b.h else port_a.h,
      if port_b.m_flow > 0 then port_b.Xi else port_a.Xi)) else 
         0 "Results in numeric Jacobian";
  */
    
  end PartialTransportIsentropicAA;
  
  redeclare replaceable partial model extends PartialTransportIsentropicAB 
    "Partial isentropic element transporting fluid between two ports without storing mass or energy (a Port_a and Port_b each, allowed in this approach)" 
    
    import Modelica_Fluid.Types;
    
    Medium.SpecificEnthalpy h_a_outflow = if flowDirection==Types.FlowDirection.Bidirectional then port_b.h - eta_ise*(port_b.h - Medium.isentropicEnthalpy(port_a.p, medium_nonDesignDirection.state)) else port_b.h;
    Medium.SpecificEnthalpy h_b_outflow = port_a.h - eta_ise*(port_a.h - Medium.isentropicEnthalpy(port_b.p, medium_designDirection.state));
    
  equation 
  /* Handle reverse and zero flow */
    port_a.H_flow = semiLinear(
            port_a.m_flow,
            port_a.h,
            h_a_outflow);
    port_b.H_flow = semiLinear(port_b.m_flow,
            port_b.h,
            h_b_outflow);
    port_a.mXi_flow = semiLinear(
            port_a.m_flow,
            port_a.Xi,
            port_b.Xi);
    
  /* Mass, energy, substance mass balance */
    port_a.m_flow + port_b.m_flow = 0;
    port_a.H_flow + port_b.H_flow + P_mechanical = 0;
    port_a.mXi_flow + port_b.mXi_flow = zeros(Medium.nXi);
    
  // Design direction of mass flow rate
    m_flow = port_a.m_flow;
    
  // Pressure difference between ports
    dp = port_a.p - port_b.p;
    
    // This approach does not provide both potential upstream properties 
    // (only upstream properties for current mass flow direction and downstream mixing properties)
    p_designDirection = port_a.p 
      "Upstream pressure if flow was in design direction";
    h_designDirection = port_a.h 
      "Approximation of upstream specific enthalpy if flow was in design direction";
    Xi_designDirection = port_a.Xi 
      "Approxiamtion of upstream mass fractions if flow was in design direction";
    p_nonDesignDirection = port_b.p 
      "Upstream pressure if flow was in non-design direction";
    h_nonDesignDirection = port_b.h 
      "Approximation of upstream specific enthalpy if flow was in non-design direction";
    Xi_nonDesignDirection = port_b.Xi 
      "Approximation of upstream mass fractions if flow was in non-design direction";
    
    // sensors
    calc_T_a = if provide_T_a then calc_T_a_medium.T else 0;
    calc_T_b = if provide_T_b then calc_T_b_medium.T else 0;
    calc_p_a = if provide_p_a then port_a.p else 0;
    calc_p_b = if provide_p_b then port_b.p else 0;
    calc_m_flow_ab = if provide_m_flow_ab then m_flow else 0;
    
    calc_T_a_medium.p = if provide_T_a then port_a.p else Medium.p_default;
    calc_T_a_medium.h = if provide_T_a then (if port_a.m_flow > 0 then port_a.h else h_a_outflow) else Medium.h_default;
    calc_T_a_medium.Xi = if provide_T_a then (if port_a.m_flow > 0 then port_a.Xi else port_b.Xi) else zeros(Medium.nXi);
    calc_T_b_medium.p = if provide_T_b then port_b.p else Medium.p_default;
    calc_T_b_medium.h = if provide_T_b then (if port_b.m_flow > 0 then port_b.h else h_b_outflow) else Medium.h_default;
    calc_T_b_medium.Xi = if provide_T_b then (if port_b.m_flow > 0 then port_b.Xi else port_a.Xi) else zeros(Medium.nXi);
    /*
  // Simpler approach for temperature like that below results in numeric Jacobian
  calc_T_a = if provide_T_a then 
    Medium.temperature(Medium.setState_phX(
      port_a.p, if port_a.m_flow > 0 then port_a.h else port_b.h,
      if port_a.m_flow > 0 then port_a.Xi else port_b.Xi)) else 
         0 "Results in numeric Jacobian";
  calc_T_b = if provide_T_b then 
    Medium.temperature(Medium.setState_phX(
      port_b.p, if port_b.m_flow > 0 then port_b.h else port_a.h,
      if port_b.m_flow > 0 then port_b.Xi else port_a.Xi)) else 
         0 "Results in numeric Jacobian";
  */
    
  end PartialTransportIsentropicAB;
  
  redeclare replaceable partial model extends PartialIdealJunction 
    "Partial infinitesimal junction model" 
    
  equation 
    connect(port_1, port_3) annotation (points=[-100,0; 0,0; 0,100], style(
      color=69, rgbcolor={0,127,255}));
    connect(port_1, port_2) annotation (points=[-100,0; 100,0], style(color=
            69,
      rgbcolor={0,127,255}));
  end PartialIdealJunction;
  
  redeclare replaceable partial model extends PartialIdealJunctionAAB 
    "Partial infinitesimal junction model (two PortA's, one PortB, not supported for all interfaces)" 
    
  equation 
    connect(port_1, port_3) annotation (points=[-100,0; 0,0; 0,100], style(
      color=69, rgbcolor={0,127,255}));
    connect(port_1, port_2) annotation (points=[-100,0; 100,0], style(color=
            69,
      rgbcolor={0,127,255}));
  end PartialIdealJunctionAAB;
  
  redeclare replaceable partial model extends PartialSource_A 
    "Partial source model with a Port_a" 
    
  equation 
    port.p = medium.p;
    port.H_flow = semiLinear(
            port.m_flow,
            port.h,
            medium.h);
    port.mXi_flow = semiLinear(
            port.m_flow,
            port.Xi,
            medium.Xi);
    
  // Interface to generic implementation
    port.m_flow = port_m_flow;
  end PartialSource_A;
  
  redeclare replaceable partial model extends PartialSource_B 
    "Partial source model with a Port_b" 
    
  equation 
    port.p = medium.p;
    port.H_flow = semiLinear(
            port.m_flow,
            port.h,
            medium.h);
    port.mXi_flow = semiLinear(
            port.m_flow,
            port.Xi,
            medium.Xi);
    
  // Interface to generic implementation
    port.m_flow = port_m_flow;
  end PartialSource_B;
  
  redeclare replaceable partial model extends PartialAsymmetricDistributedPipe 
  // Taken one-to-one from Modelica_Fluid
    
  // Covered herein:
  //   Properties associated with ports
  //   Aliases to connector variables
  //   Connector equations
    
    extends FluidDiscretization.PartialDistributedFlow;
    
  equation 
  // Aliases to connectors
    port_a_p = port_a.p;
    port_b_p = port_b.p;
    port_a.m_flow = m_flow[1];
    port_b.m_flow = -m_flow[n + 1];
    H_flow[1] = port_a.H_flow;
    H_flow[n + 1] = -port_b.H_flow;
    mXi_flow[1, :] = port_a.mXi_flow;
    mXi_flow[n + 1, :] = -port_b.mXi_flow;
  // Quantities related to connectors
    v[1] = m_flow[1]/d_a/area;
    v[n + 1] = m_flow[n + 1]/d_b/area;
    // For pressure drop correlations
    eta_a = 0 "This is only used for the symmetric discretizations";
    d_a = 0 "This is only used for the symmetric discretizations";
    eta_b = if not WallFriction.use_eta then 1.e-10 else (if 
      use_eta_nominal then eta_nominal else (if use_approxPortProperties then 
            eta[n] else (if m_flow[n + 1] < 0 then Medium.dynamicViscosity(
      Medium.setState_phX(
            port_b.p,
            port_b.h,
            port_b.Xi)) else eta[n])));
    d_b = if use_d_nominal then d_nominal else (if use_approxPortProperties then 
            d[n] else (if m_flow[n + 1] >= 0 then d[n] else 
      Medium.density_phX(
            port_b.p,
            port_b.h,
            port_b.Xi)));
    
  // Connector equations
  // Enthalpy flow
    port_a.H_flow = semiLinear(
            port_a.m_flow,
            port_a.h,
            medium[1].h);
    port_b.H_flow = semiLinear(
            port_b.m_flow,
            port_b.h,
            medium[n].h);
  // Substance mass flow
    port_a.mXi_flow = semiLinear(
            port_a.m_flow,
            port_a.Xi,
            medium[1].Xi);
    port_b.mXi_flow = semiLinear(
            port_b.m_flow,
            port_b.Xi,
            medium[n].Xi);
    
    annotation (
      Diagram,
      Icon,
      Documentation(info="<html>
<p>The model <b>PartialDistributedFlow</b> is used as a base class for pipe flows with one-dimensional spatial discretization according to the finite volume method. The flow path is divided into <tt><b>n</b></tt> segments.</p>
 
<p><b>Mass and energy balances</b></p>
<p>One total mass and one energy balance is formed across each segment. If the medium contains more than one component, substance mass balances are added. Changes in potential and kinetic energy are neglected in the energy balance. The following source (or sink) terms are used in the balances and must be specified in extending models to complete this partial class:</p>
<ul>
<li>Energy balance: <tt><b>Qs_flow</b></tt>, e.g. convective or latent heat flow rate across segment boundary, and <tt><b>Ws_flow</b></tt>, e.g. mechanical power</li>
<li>Total mass balance: <tt><b>ms_flow</b></tt>, e.g. condensing mass flow of negligible volume such as water in moist air</li>
<li>Substance mass balance: <tt><b>msXi_flow</b></tt>, as above</li>
</ul>
If the flag <tt>static</tt> is <b>true</b> then no mass or energy is stored in the component and the mass and energy balances are reduced to a quasi steady-state formulation. It should be noted that dynamic balances are required if flow reversal should be allowed.
<p>In addition the volume vector <tt><b>Vi</b></tt>, which specifies the volume of each segment and the pressure drop (or rise) in each segment <tt><b>dp</b></tt> must be provided in the extending class.
 
<p><b>Momentum balance</b></p>
<p>The momentum balance is always static, i.e. no dynamic momentum term is used. The momentum balances are formed across the segment boundaries (staggered grid). The default symmetric model is characterized by half a momentum balance on each end of the flow model resulting in a total of n-1 full and 2 half momentum balances. Connecting two pipes therefore results in an algebraic pressure at the ports. Specifying a good start value for the port pressure is essential in order to solve large systems. Non-symmetric variations are obtained by chosing a different value for the parameter <tt><b>modelStructure</b></tt>. Options include:
<ul>
<li><tt>a_v_b</tt>: default setting with two half momentum balances</li>
<li><tt>av_b</tt>: full momentum balance between nth volume and <tt>port_b</tt>, potential pressure state at <tt>port_a</tt></li>
<li><tt>a_vb</tt>: full momentum balance between first volume and <tt>port_a</tt>, potential pressure state at <tt>port_b</tt></li>
<li><tt>avb</tt>: n-1 momentum balances between first and nth volume, potential pressure states at both ports. It's use should be avoided, since not the entire pipe length is taken into account.
</ul></p>
 
<p>The term <tt>dp</tt> is unspecified in this partial class. When extending from this model it may contain
<ul>
<li>pressure drop due to friction and other dissipative losses</li>
<li>changes in pressure resulting from significant variation of flow velocity along the flow path (with the assumption of a constant cross sectional area it must result from fluid density changes, such as in two-phase flow)</li>
<li>gravity effects for non-horizontal pipes</li>
</ul>
At least one relationship between pressure difference and massflow rate (dp=dp(m_flow)) is required in the extending class for this term to receive a fully determined model.
 
When connecting two components, e.g. two pipes, the momentum balance across the connection point reduces to</p> 
<pre>pipe1.port_b.p = pipe2.port_a.p</pre>
<p>This is only true if the flow velocity remains the same on each side of the connection. For any significant change in diameter (and if the resulting effects, such as change in kinetic energy, cannot be neglected) an adapter component should be used. This also allows for taking into account friction losses with respect to the actual geometry of the connection point.</p>
 
</html>",
      revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
    
  end PartialAsymmetricDistributedPipe;
  
  redeclare replaceable partial model extends PartialSymmetricDistributedPipe 
  // Taken one-to-one from Modelica_Fluid
    
  // Covered herein:
  //   Properties associated with ports
  //   Aliases to connector variables
  //   Connector equations
    
    extends FluidDiscretization.PartialDistributedFlow;
    
  equation 
  // Aliases to connectors
    port_a_p = port_a.p;
    port_b_p = port_b.p;
    port_a.m_flow = m_flow[1];
    port_b.m_flow = -m_flow[n + 1];
    H_flow[1] = port_a.H_flow;
    H_flow[n + 1] = -port_b.H_flow;
    mXi_flow[1, :] = port_a.mXi_flow;
    mXi_flow[n + 1, :] = -port_b.mXi_flow;
  // Quantities related to connectors
    v[1] = m_flow[1]/d_a/area;
    v[n + 1] = m_flow[n + 1]/d_b/area;
    eta_a = if not WallFriction.use_eta then 1.e-10 else (if 
      use_eta_nominal then eta_nominal else (if use_approxPortProperties then 
            eta[1] else (if m_flow[1] >= 0 then Medium.dynamicViscosity(
      Medium.setState_phX(
            port_a.p,
            port_a.h,
            port_a.Xi)) else eta[1])));
    eta_b = if not WallFriction.use_eta then 1.e-10 else (if 
      use_eta_nominal then eta_nominal else (if use_approxPortProperties then 
            eta[n] else (if m_flow[n + 1] < 0 then Medium.dynamicViscosity(
      Medium.setState_phX(
            port_b.p,
            port_b.h,
            port_b.Xi)) else eta[n])));
    d_a = if use_d_nominal then d_nominal else (if use_approxPortProperties then 
            d[1] else (if m_flow[1] >= 0 then Medium.density_phX(
            port_a.p,
            port_a.h,
            port_a.Xi) else d[1]));
    d_b = if use_d_nominal then d_nominal else (if use_approxPortProperties then 
            d[n] else (if m_flow[n + 1] >= 0 then d[n] else 
      Medium.density_phX(
            port_b.p,
            port_b.h,
            port_b.Xi)));
    
  // Connector equations
  // Enthalpy flow
    port_a.H_flow = semiLinear(
            port_a.m_flow,
            port_a.h,
            medium[1].h);
    port_b.H_flow = semiLinear(
            port_b.m_flow,
            port_b.h,
            medium[n].h);
  // Substance mass flow
    port_a.mXi_flow = semiLinear(
            port_a.m_flow,
            port_a.Xi,
            medium[1].Xi);
    port_b.mXi_flow = semiLinear(
            port_b.m_flow,
            port_b.Xi,
            medium[n].Xi);
    
    annotation (
      Diagram,
      Icon,
      Documentation(info="<html>
<p>The model <b>PartialDistributedFlow</b> is used as a base class for pipe flows with one-dimensional spatial discretization according to the finite volume method. The flow path is divided into <tt><b>n</b></tt> segments.</p>
 
<p><b>Mass and energy balances</b></p>
<p>One total mass and one energy balance is formed across each segment. If the medium contains more than one component, substance mass balances are added. Changes in potential and kinetic energy are neglected in the energy balance. The following source (or sink) terms are used in the balances and must be specified in extending models to complete this partial class:</p>
<ul>
<li>Energy balance: <tt><b>Qs_flow</b></tt>, e.g. convective or latent heat flow rate across segment boundary, and <tt><b>Ws_flow</b></tt>, e.g. mechanical power</li>
<li>Total mass balance: <tt><b>ms_flow</b></tt>, e.g. condensing mass flow of negligible volume such as water in moist air</li>
<li>Substance mass balance: <tt><b>msXi_flow</b></tt>, as above</li>
</ul>
If the flag <tt>static</tt> is <b>true</b> then no mass or energy is stored in the component and the mass and energy balances are reduced to a quasi steady-state formulation. It should be noted that dynamic balances are required if flow reversal should be allowed.
<p>In addition the volume vector <tt><b>Vi</b></tt>, which specifies the volume of each segment and the pressure drop (or rise) in each segment <tt><b>dp</b></tt> must be provided in the extending class.
 
<p><b>Momentum balance</b></p>
<p>The momentum balance is always static, i.e. no dynamic momentum term is used. The momentum balances are formed across the segment boundaries (staggered grid). The default symmetric model is characterized by half a momentum balance on each end of the flow model resulting in a total of n-1 full and 2 half momentum balances. Connecting two pipes therefore results in an algebraic pressure at the ports. Specifying a good start value for the port pressure is essential in order to solve large systems. Non-symmetric variations are obtained by chosing a different value for the parameter <tt><b>modelStructure</b></tt>. Options include:
<ul>
<li><tt>a_v_b</tt>: default setting with two half momentum balances</li>
<li><tt>av_b</tt>: full momentum balance between nth volume and <tt>port_b</tt>, potential pressure state at <tt>port_a</tt></li>
<li><tt>a_vb</tt>: full momentum balance between first volume and <tt>port_a</tt>, potential pressure state at <tt>port_b</tt></li>
<li><tt>avb</tt>: n-1 momentum balances between first and nth volume, potential pressure states at both ports. It's use should be avoided, since not the entire pipe length is taken into account.
</ul></p>
 
<p>The term <tt>dp</tt> is unspecified in this partial class. When extending from this model it may contain
<ul>
<li>pressure drop due to friction and other dissipative losses</li>
<li>changes in pressure resulting from significant variation of flow velocity along the flow path (with the assumption of a constant cross sectional area it must result from fluid density changes, such as in two-phase flow)</li>
<li>gravity effects for non-horizontal pipes</li>
</ul>
At least one relationship between pressure difference and massflow rate (dp=dp(m_flow)) is required in the extending class for this term to receive a fully determined model.
 
When connecting two components, e.g. two pipes, the momentum balance across the connection point reduces to</p> 
<pre>pipe1.port_b.p = pipe2.port_a.p</pre>
<p>This is only true if the flow velocity remains the same on each side of the connection. For any significant change in diameter (and if the resulting effects, such as change in kinetic energy, cannot be neglected) an adapter component should be used. This also allows for taking into account friction losses with respect to the actual geometry of the connection point.</p>
 
</html>",
      revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  end PartialSymmetricDistributedPipe;
  
end SemiLinear_A1;
