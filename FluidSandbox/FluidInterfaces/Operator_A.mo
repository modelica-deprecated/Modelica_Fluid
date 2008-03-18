within FluidSandbox.FluidInterfaces;
package Operator_A 
  "Implementation A using a special operator inflow() (available in Dymola 7.0 beta 3x)" 
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
    
    // This approach provides both potential upstream properties (independent of current mass flow direction)
    p_designDirection = port_a.p 
      "Upstream pressure if flow was in design direction";
    h_designDirection = inflow(port_a.m_flow, port_a.h) 
      "Upstream specific enthalpy if flow was in design direction";
    Xi_designDirection =                           port_a.Xi 
      "Upstream mass fractions if flow was in design direction";
                         /*inflow(port_a.m_flow, */
                                                            //) 
    p_nonDesignDirection = port_b.p 
      "Upstream pressure if flow was in non-design direction";
    h_nonDesignDirection = inflow(port_b.m_flow, port_b.h) 
      "Upstream specific enthalpy if flow was in non-design direction";
    Xi_nonDesignDirection =                           port_b.Xi 
      "Upstream mass fractions if flow was in non-design direction";
                            /*inflow(port_b.m_flow, */
                                                               //) 
    
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
    
    // This approach provides both potential upstream properties (independent of current mass flow direction)
    p_designDirection = port_a.p 
      "Upstream pressure if flow was in design direction";
    h_designDirection = inflow(port_a.m_flow, port_a.h) 
      "Upstream specific enthalpy if flow was in design direction";
    Xi_designDirection =                           port_a.Xi 
      "Upstream mass fractions if flow was in design direction";
                         /*inflow(port_a.m_flow, */
                                                            //) 
    p_nonDesignDirection = port_b.p 
      "Upstream pressure if flow was in non-design direction";
    h_nonDesignDirection = inflow(port_b.m_flow, port_b.h) 
      "Upstream specific enthalpy if flow was in non-design direction";
    Xi_nonDesignDirection =                           port_b.Xi 
      "Upstream mass fractions if flow was in non-design direction";
                            /*inflow(port_b.m_flow, */
                                                               //) 
    
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
    
    // This approach provides both potential upstream properties (independent of current mass flow direction)
    p_designDirection = port_a.p 
      "Upstream pressure if flow was in design direction";
    h_designDirection = inflow(port_a.m_flow, port_a.h) 
      "Upstream specific enthalpy if flow was in design direction";
    Xi_designDirection =                           port_a.Xi 
      "Upstream mass fractions if flow was in design direction";
                         /*inflow(port_a.m_flow, */
                                                            //) 
    p_nonDesignDirection = port_b.p 
      "Upstream pressure if flow was in non-design direction";
    h_nonDesignDirection = inflow(port_b.m_flow, port_b.h) 
      "Upstream specific enthalpy if flow was in non-design direction";
    Xi_nonDesignDirection =                           port_b.Xi 
      "Upstream mass fractions if flow was in non-design direction";
                            /*inflow(port_b.m_flow, */
                                                               //) 
    
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
    
  end PartialTransportIsenthalpicAB;
  
  redeclare replaceable partial model extends PartialTransportIsentropic 
    "Partial isentropic element transporting fluid between two ports without storing mass or energy (two Port_b's)" 
    
    import Modelica_Fluid.Types;
    
    Medium.SpecificEnthalpy h_a_outflow = if flowDirection==Types.FlowDirection.Bidirectional then h_nonDesignDirection - eta_ise*(h_nonDesignDirection - Medium.isentropicEnthalpy(port_a.p, medium_nonDesignDirection.state)) else h_nonDesignDirection;
    Medium.SpecificEnthalpy h_b_outflow = h_designDirection - eta_ise*(h_designDirection - Medium.isentropicEnthalpy(port_b.p, medium_designDirection.state));
    
    /*
  // Isentropic process
  // Preferably, if supported by Medium model
  h_ba_isentropic = Medium.isentropicEnthalpy(port_a.p, medium_nonDesignDirection.state);
  h_ab_isentropic = Medium.isentropicEnthalpy(port_b.p, medium_designDirection.state);
  // Implicit equation if not supported
  Medium.specificEntropy(medium_nonDesignDirection) = Medium.specificEntropy(Medium.setState_phX(port_a.p, h_ba_isentropic, port_b.Xi));
  Medium.specificEntropy(medium_designDirection) = Medium.specificEntropy(Medium.setState_phX(port_b.p, h_ab_isentropic, port_a.Xi));
  */
    
  equation 
    /* Handle reverse and zero flow */
    port_a.H_flow = semiLinear(
      port_a.m_flow,
      port_a.h,
      h_a_outflow) 
      "According to Sven Erik, semiLinear has to be provided a variable, not an expression";
    port_b.H_flow = semiLinear(
      port_b.m_flow,
      port_b.h,
      h_b_outflow) 
      "According to Sven Erik, semiLinear has to be provided a variable, not an expression";
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
    
    // This approach provides both potential upstream properties (independent of current mass flow direction)
    p_designDirection = port_a.p 
      "Upstream pressure if flow was in design direction";
    h_designDirection = inflow(port_a.m_flow, port_a.h) 
      "Upstream specific enthalpy if flow was in design direction";
    Xi_designDirection = inflow(port_a.m_flow, port_a.Xi) 
      "Upstream mass fractions if flow was in design direction";
    
    p_nonDesignDirection = port_b.p 
      "Upstream pressure if flow was in non-design direction";
    h_nonDesignDirection = inflow(port_b.m_flow, port_b.h) 
      "Upstream specific enthalpy if flow was in non-design direction";
    Xi_nonDesignDirection = inflow(port_b.m_flow, port_b.Xi) 
      "Upstream mass fractions if flow was in non-design direction";
    
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
    
  end PartialTransportIsentropic;
  
  redeclare replaceable partial model extends PartialTransportIsentropicAA 
    "Partial isentropic element transporting fluid between two ports without storing mass or energy (two Port_a's, allowed in this approach)" 
    
    import Modelica_Fluid.Types;
    
    Medium.SpecificEnthalpy h_a_outflow = if flowDirection==Types.FlowDirection.Bidirectional then h_nonDesignDirection - eta_ise*(h_nonDesignDirection - Medium.isentropicEnthalpy(port_a.p, medium_nonDesignDirection.state)) else h_nonDesignDirection;
    Medium.SpecificEnthalpy h_b_outflow = h_designDirection - eta_ise*(h_designDirection - Medium.isentropicEnthalpy(port_b.p, medium_designDirection.state));
    
    /*
  // Isentropic process
  // Preferably, if supported by Medium model
  h_ba_isentropic = Medium.isentropicEnthalpy(port_a.p, medium_nonDesignDirection.state);
  h_ab_isentropic = Medium.isentropicEnthalpy(port_b.p, medium_designDirection.state);
  // Implicit equation if not supported
  Medium.specificEntropy(medium_nonDesignDirection) = Medium.specificEntropy(Medium.setState_phX(port_a.p, h_ba_isentropic, port_b.Xi));
  Medium.specificEntropy(medium_designDirection) = Medium.specificEntropy(Medium.setState_phX(port_b.p, h_ab_isentropic, port_a.Xi));
  */
    
  equation 
    /* Handle reverse and zero flow */
    port_a.H_flow = semiLinear(
      port_a.m_flow,
      port_a.h,
      h_a_outflow) 
      "According to Sven Erik, semiLinear has to be provided a variable, not an expression";
    port_b.H_flow = semiLinear(
      port_b.m_flow,
      port_b.h,
      h_b_outflow) 
      "According to Sven Erik, semiLinear has to be provided a variable, not an expression";
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
    
    // This approach provides both potential upstream properties (independent of current mass flow direction)
    p_designDirection = port_a.p 
      "Upstream pressure if flow was in design direction";
    h_designDirection = inflow(port_a.m_flow, port_a.h) 
      "Upstream specific enthalpy if flow was in design direction";
    Xi_designDirection = inflow(port_a.m_flow, port_a.Xi) 
      "Upstream mass fractions if flow was in design direction";
    
    p_nonDesignDirection = port_b.p 
      "Upstream pressure if flow was in non-design direction";
    h_nonDesignDirection = inflow(port_b.m_flow, port_b.h) 
      "Upstream specific enthalpy if flow was in non-design direction";
    Xi_nonDesignDirection = inflow(port_b.m_flow, port_b.Xi) 
      "Upstream mass fractions if flow was in non-design direction";
    
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
    
  end PartialTransportIsentropicAA;
  
  redeclare replaceable partial model extends PartialTransportIsentropicAB 
    "Partial isentropic element transporting fluid between two ports without storing mass or energy (a Port_a and Port_b each, allowed in this approach)" 
    
    import Modelica_Fluid.Types;
    
    Medium.SpecificEnthalpy h_a_outflow = if flowDirection==Types.FlowDirection.Bidirectional then h_nonDesignDirection - eta_ise*(h_nonDesignDirection - Medium.isentropicEnthalpy(port_a.p, medium_nonDesignDirection.state)) else h_nonDesignDirection;
    Medium.SpecificEnthalpy h_b_outflow = h_designDirection - eta_ise*(h_designDirection - Medium.isentropicEnthalpy(port_b.p, medium_designDirection.state));
    
    /*
  // Isentropic process
  // Preferably, if supported by Medium model
  h_ba_isentropic = Medium.isentropicEnthalpy(port_a.p, medium_nonDesignDirection.state);
  h_ab_isentropic = Medium.isentropicEnthalpy(port_b.p, medium_designDirection.state);
  // Implicit equation if not supported
  Medium.specificEntropy(medium_nonDesignDirection) = Medium.specificEntropy(Medium.setState_phX(port_a.p, h_ba_isentropic, port_b.Xi));
  Medium.specificEntropy(medium_designDirection) = Medium.specificEntropy(Medium.setState_phX(port_b.p, h_ab_isentropic, port_a.Xi));
  */
    
  equation 
    /* Handle reverse and zero flow */
    port_a.H_flow = semiLinear(
      port_a.m_flow,
      port_a.h,
      h_a_outflow) 
      "According to Sven Erik, semiLinear has to be provided a variable, not an expression";
    port_b.H_flow = semiLinear(
      port_b.m_flow,
      port_b.h,
      h_b_outflow) 
      "According to Sven Erik, semiLinear has to be provided a variable, not an expression";
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
    
    // This approach provides both potential upstream properties (independent of current mass flow direction)
    p_designDirection = port_a.p 
      "Upstream pressure if flow was in design direction";
    h_designDirection = inflow(port_a.m_flow, port_a.h) 
      "Upstream specific enthalpy if flow was in design direction";
    Xi_designDirection = inflow(port_a.m_flow, port_a.Xi) 
      "Upstream mass fractions if flow was in design direction";
    
    p_nonDesignDirection = port_b.p 
      "Upstream pressure if flow was in non-design direction";
    h_nonDesignDirection = inflow(port_b.m_flow, port_b.h) 
      "Upstream specific enthalpy if flow was in non-design direction";
    Xi_nonDesignDirection = inflow(port_b.m_flow, port_b.Xi) 
      "Upstream mass fractions if flow was in non-design direction";
    
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
    
  end PartialTransportIsentropicAB;
  
  redeclare replaceable partial model extends PartialIdealJunction 
    "Partial infinitesimal junction model" 
    
  equation 
    connect(port_1, port_3) 
                          annotation (points=[-100,0; 0,0; 0,100], style(
    color=69, rgbcolor={0,127,255}));
    connect(port_1, port_2) 
                          annotation (points=[-100,0; 100,0], style(color=
          69,
    rgbcolor={0,127,255}));
  end PartialIdealJunction;
  
  redeclare replaceable partial model extends PartialIdealJunctionAAB 
    "Partial infinitesimal junction model (two PortA's, one PortB, not supported for all interfaces)" 
    
  equation 
    connect(port_1, port_3) 
                          annotation (points=[-100,0; 0,0; 0,100], style(
    color=69, rgbcolor={0,127,255}));
    connect(port_1, port_2) 
                          annotation (points=[-100,0; 100,0], style(color=
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
end Operator_A;
