package EffortPairs_A 
  "Implementation A using pairs of effort variables on connectors (similar to the ThermoPower library)" 
  extends Interfaces.PartialFluidInterface(usesNewConnectionSemantics=false);
  redeclare replaceable connector extends FluidPort 
    "Interface for quasi one-dimensional fluid flow in a piping network (incompressible or compressible, one or more phases, one or more substances)" 
    
    Medium.AbsolutePressure p "Pressure in the connection point";
    flow Medium.MassFlowRate m_flow 
      "Mass flow rate from the connection point into the component";
    
  end FluidPort;
  
  redeclare replaceable connector extends FluidPort_a 
    "Generic fluid connector at design inlet" 
    
    output Medium.SpecificEnthalpy h_a 
      "Specific enthalpy in PortA (upstream if mass flows from A to B, downstream otherwise)";
    input Medium.SpecificEnthalpy h_b 
      "Specific enthalpy in PortB (upstream if mass flows from B to A, downstream otherwise)";
    
    output Medium.MassFraction Xi_a[Medium.nXi] 
      "Independent substance mass fractions m_i/m in PortA (upstream if mass flows from A to B, downstream otherwise)";
    input Medium.MassFraction Xi_b[Medium.nXi] 
      "Independent substance mass fractions m_i/m in PortB (upstream if mass flows from B to A, downstream otherwise)";
    
  end FluidPort_a;
  
  redeclare replaceable connector extends FluidPort_b 
    "Generic fluid connector at design outlet" 
    
    input Medium.SpecificEnthalpy h_a 
      "Specific enthalpy in PortA (upstream if mass flows from A to B, downstream otherwise)";
    output Medium.SpecificEnthalpy h_b 
      "Specific enthalpy in PortB (upstream if mass flows from B to A, downstream otherwise)";
    
    input Medium.MassFraction Xi_a[Medium.nXi] 
      "Independent substance mass fractions m_i/m in PortA (upstream if mass flows from A to B, downstream otherwise)";
    output Medium.MassFraction Xi_b[Medium.nXi] 
      "Independent substance mass fractions m_i/m in PortB (upstream if mass flows from B to A, downstream otherwise)";
    
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
    port_a.h_a = medium.h*ones(n_a);
    port_b.h_a = medium.h*ones(n_b) "Port_a instance, too";
    
    // Substance mass flows
    for j in 1:Medium.nXi loop
      port_a[:].Xi_a[j] = medium.Xi[j]*ones(n_a);
      port_b[:].Xi_a[j] = medium.Xi[j]*ones(n_b) 
        "port_b is a Port_a instance, too";
    end for;
    
    // Net flow rates
    m_flow_net = sum(port_a.m_flow) + sum(port_b.m_flow);
    for j in 1:Medium.nXi loop
      mXi_flow_net[j] = sum({port_a[i].m_flow*(if port_a[i].m_flow > 0 then 
              port_a[i].Xi_b[j] else port_a[i].Xi_a[j]) for i in 1:n_a}) +
        sum({port_b[i].m_flow*(if port_b[i].m_flow > 0 then port_b[i].Xi_b[
        j] else port_b[i].Xi_a[j]) for i in 1:n_b});
    end for;
    H_flow_net = sum({port_a[i].m_flow*(if port_a[i].m_flow > 0 then port_a[
      i].h_b else port_a[i].h_a) for i in 1:n_a}) + sum({port_b[i].m_flow*(
      if port_b[i].m_flow > 0 then port_b[i].h_b else port_b[i].h_a) for i in 
          1:n_b});
  /*
  All substances at once
  mXi_flow_net = {sum(port_a[:].mXi_flow[i]) + sum(port_b[:].mXi_flow[i]) for i in 1:Medium.nXi};
 
*/
    
  end PartialLumpedVolume;
  
  redeclare replaceable partial model extends PartialTransportIsenthalpic 
    "Partial isenthalpic element transporting fluid between two ports without storing mass or energy (two Port_b's)" 
    
  equation 
    // Mass balance
    port_a.m_flow + port_b.m_flow = 0;
    // Enthalpy propagation, energy balance
    port_a.h_b = port_b.h_a;
    port_b.h_b = port_a.h_a;
    // Mass fraction propagation, substance mass balance
    port_a.Xi_b = port_b.Xi_a;
    port_b.Xi_b = port_a.Xi_a;
    
    // Design direction of mass flow rate
    m_flow = port_a.m_flow;
    
    // Pressure difference between ports
    dp = port_a.p - port_b.p;
    
    // This approach provides upstream and downstream properties
    p_designDirection = port_a.p 
      "Upstream pressure if flow is in design direction";
    h_designDirection = port_a.h_a 
      "Upstream specific enthalpy if flow is in design direction";
    Xi_designDirection = port_a.Xi_a 
      "Upstream mass fractions if flow is in design direction";
    p_nonDesignDirection = port_b.p 
      "Upstream pressure if flow is in non-design direction";
    h_nonDesignDirection = port_b.h_a 
      "Upstream specific enthalpy if flow is in non-design direction";
    Xi_nonDesignDirection = port_b.Xi_a 
      "Upstream mass fractions if flow is in non-design direction";
    
  end PartialTransportIsenthalpic;
  
  redeclare replaceable partial model extends PartialTransportIsenthalpicAA 
    "Partial isenthalpic element transporting fluid between two ports without storing mass or energy (two Port_a's, allowed in this approach)" 
    
  equation 
    // Mass balance
    port_a.m_flow + port_b.m_flow = 0;
    // Enthalpy propagation, energy balance
    port_a.h_a = port_b.h_b;
    port_b.h_a = port_a.h_b;
    // Mass fraction propagation, substance mass balance
    port_a.Xi_a = port_b.Xi_b;
    port_b.Xi_a = port_a.Xi_b;
    
    // Design direction of mass flow rate
    m_flow = port_a.m_flow;
    
    // Pressure difference between ports
    dp = port_a.p - port_b.p;
    
    // This approach provides upstream and downstream properties
    p_designDirection = port_a.p 
      "Upstream pressure if flow is in design direction";
    h_designDirection = port_a.h_b 
      "Upstream specific enthalpy if flow is in design direction";
    Xi_designDirection = port_a.Xi_b 
      "Upstream mass fractions if flow is in design direction";
    p_nonDesignDirection = port_b.p 
      "Upstream pressure if flow is in non-design direction";
    h_nonDesignDirection = port_b.h_b 
      "Upstream specific enthalpy if flow is in non-design direction";
    Xi_nonDesignDirection = port_b.Xi_b 
      "Upstream mass fractions if flow is in non-design direction";
    
  end PartialTransportIsenthalpicAA;
  
  redeclare replaceable partial model extends PartialTransportIsenthalpicAB 
    "Partial isenthalpic element transporting fluid between two ports without storing mass or energy (a Port_a and Port_b each, allowed in this approach)" 
    
  equation 
    // Mass balance
    port_a.m_flow + port_b.m_flow = 0;
    // Enthalpy propagation, energy balance
    port_a.h_a = port_b.h_a;
    port_b.h_b = port_a.h_b;
    // Mass fraction propagation, substance mass balance
    port_a.Xi_a = port_b.Xi_a;
    port_b.Xi_b = port_a.Xi_b;
    
    // Design direction of mass flow rate
    m_flow = port_a.m_flow;
    
    // Pressure difference between ports
    dp = port_a.p - port_b.p;
    
    // This approach provides upstream and downstream properties
    p_designDirection = port_a.p 
      "Upstream pressure if flow is in design direction";
    h_designDirection = port_a.h_b 
      "Upstream specific enthalpy if flow is in design direction";
    Xi_designDirection = port_a.Xi_b 
      "Upstream mass fractions if flow is in design direction";
    p_nonDesignDirection = port_b.p 
      "Upstream pressure if flow is in non-design direction";
    h_nonDesignDirection = port_b.h_a 
      "Upstream specific enthalpy if flow is in non-design direction";
    Xi_nonDesignDirection = port_b.Xi_a 
      "Upstream mass fractions if flow is in non-design direction";
    
  end PartialTransportIsenthalpicAB;
  
  redeclare replaceable partial model extends PartialTransportIsentropic 
    "Partial isentropic element transporting fluid between two ports without storing mass or energy (two Port_b's)" 
    
  equation 
    // Mass balance
    port_a.m_flow + port_b.m_flow = 0;
    
    // Enthalpy propagation, energy balance
    port_b.h_b = port_a.h_a - eta_ise*(port_a.h_a - Medium.isentropicEnthalpy(port_b.p, medium_a)) 
      "Design mass flow direction";
    port_a.h_b = port_b.h_a - eta_ise*(port_b.h_a - Medium.isentropicEnthalpy(port_a.p, medium_b)) 
      "Non-design mass flow direction";
    
    P_mechanical = -noEvent((port_a.m_flow*(if port_a.m_flow > 0 then port_a.h_a else port_a.h_b)) + (port_b.m_flow*(
      if port_b.m_flow > 0 then port_b.h_a else port_b.h_b)));
    
    // Mass fraction propagation, substance mass balance
    port_a.Xi_b = port_b.Xi_a;
    port_b.Xi_b = port_a.Xi_a;
    
    // Design direction of mass flow rate
    m_flow = port_a.m_flow;
    
    // Pressure difference between ports
    dp = port_a.p - port_b.p;
    
    // This approach provides upstream and downstream properties
    p_designDirection = port_a.p 
      "Upstream pressure if flow is in design direction";
    h_designDirection = port_a.h_a 
      "Upstream specific enthalpy if flow is in design direction";
    Xi_designDirection = port_a.Xi_a 
      "Upstream mass fractions if flow is in design direction";
    p_nonDesignDirection = port_b.p 
      "Upstream pressure if flow is in non-design direction";
    h_nonDesignDirection = port_b.h_a 
      "Upstream specific enthalpy if flow is in non-design direction";
    Xi_nonDesignDirection = port_b.Xi_a 
      "Upstream mass fractions if flow is in non-design direction";
    
  end PartialTransportIsentropic;

  redeclare replaceable partial model extends PartialTransportIsentropicAA 
    "Partial isentropic element transporting fluid between two ports without storing mass or energy (two Port_a's, allowed in this approach)" 
    
  equation 
    // Mass balance
    port_a.m_flow + port_b.m_flow = 0;
    
    // Enthalpy propagation, energy balance
    port_b.h_a = port_a.h_b - eta_ise*(port_a.h_b - Medium.isentropicEnthalpy(port_b.p, medium_a)) 
      "Design mass flow direction";
    port_a.h_a = port_b.h_b - eta_ise*(port_b.h_b - Medium.isentropicEnthalpy(port_a.p, medium_b)) 
      "Non-design mass flow direction";
    
    P_mechanical = -noEvent((port_a.m_flow*(if port_a.m_flow > 0 then port_a.h_b else port_a.h_a)) + (port_b.m_flow*(
      if port_b.m_flow > 0 then port_b.h_b else port_b.h_a)));
    
    // Mass fraction propagation, substance mass balance
    port_a.Xi_b = port_b.Xi_a;
    port_b.Xi_b = port_a.Xi_a;
    
    // Design direction of mass flow rate
    m_flow = port_a.m_flow;
    
    // Pressure difference between ports
    dp = port_a.p - port_b.p;
    
    // This approach provides upstream and downstream properties
    p_designDirection = port_a.p 
      "Upstream pressure if flow is in design direction";
    h_designDirection = port_a.h_b 
      "Upstream specific enthalpy if flow is in design direction";
    Xi_designDirection = port_a.Xi_b 
      "Upstream mass fractions if flow is in design direction";
    p_nonDesignDirection = port_b.p 
      "Upstream pressure if flow is in non-design direction";
    h_nonDesignDirection = port_b.h_b 
      "Upstream specific enthalpy if flow is in non-design direction";
    Xi_nonDesignDirection = port_b.Xi_b 
      "Upstream mass fractions if flow is in non-design direction";
    
  end PartialTransportIsentropicAA;

  redeclare replaceable partial model extends PartialTransportIsentropicAB 
    "Partial isentropic element transporting fluid between two ports without storing mass or energy (a Port_a and Port_b each, allowed in this approach)" 
    
  equation 
    // Mass balance
    port_a.m_flow + port_b.m_flow = 0;
    
    // Enthalpy propagation, energy balance
    port_b.h_b = port_a.h_b - eta_ise*(port_a.h_b - Medium.isentropicEnthalpy(port_b.p, medium_a)) 
      "Design mass flow direction";
    port_a.h_a = port_b.h_a - eta_ise*(port_b.h_a - Medium.isentropicEnthalpy(port_a.p, medium_b)) 
      "Non-design mass flow direction";
    
    P_mechanical = -noEvent((port_a.m_flow*(if port_a.m_flow > 0 then port_a.h_b else port_a.h_a)) + (port_b.m_flow*(
      if port_b.m_flow > 0 then port_b.h_a else port_b.h_b)));
    
    // Mass fraction propagation, substance mass balance
    port_a.Xi_b = port_b.Xi_a;
    port_b.Xi_b = port_a.Xi_a;
    
    // Design direction of mass flow rate
    m_flow = port_a.m_flow;
    
    // Pressure difference between ports
    dp = port_a.p - port_b.p;
    
    // This approach provides upstream and downstream properties
    p_designDirection = port_a.p 
      "Upstream pressure if flow is in design direction";
    h_designDirection = port_a.h_b 
      "Upstream specific enthalpy if flow is in design direction";
    Xi_designDirection = port_a.Xi_b 
      "Upstream mass fractions if flow is in design direction";
    p_nonDesignDirection = port_b.p 
      "Upstream pressure if flow is in non-design direction";
    h_nonDesignDirection = port_b.h_a 
      "Upstream specific enthalpy if flow is in non-design direction";
    Xi_nonDesignDirection = port_b.Xi_a 
      "Upstream mass fractions if flow is in non-design direction";
    
  end PartialTransportIsentropicAB;

  redeclare replaceable partial model extends PartialIdealJunction 
    "Partial infinitesimal junction model" 
    
    Medium.AbsolutePressure p "Pressure";
    
  equation 
    // Mass balance
    port_1.m_flow + port_2.m_flow + port_3.m_flow = 0;
    
    // Pressure
    port_1.p = p;
    port_2.p = p;
    port_3.p = p;
    
    // Flow directions
    if noEvent(port_2.m_flow > 0 and port_3.m_flow < 0) then
      port_1.h_a = port_2.h_b;
      port_1.Xi_a = port_2.Xi_b;
    elseif noEvent(port_2.m_flow < 0 and port_3.m_flow > 0) then
      port_1.h_a = port_3.h_b;
      port_1.Xi_a = port_3.Xi_b;
    elseif noEvent(port_2.m_flow > 0 and port_3.m_flow > 0) then
      port_1.h_a = (port_2.m_flow*port_2.h_b + port_3.m_flow*port_3.h_b)/(port_2.m_flow
         + port_3.m_flow);
      port_1.Xi_a = (port_2.m_flow*port_2.Xi_b + port_3.m_flow*port_3.Xi_b)/(port_2.m_flow
         + port_3.m_flow);
    else
      port_1.h_a = port_1.h_b 
        "The theoretical inverse flow direction means in this case that all three flows exit the junction";
      port_1.Xi_a = port_1.Xi_b 
        "The theoretical inverse flow direction means in this case that all three flows exit the junction";
    end if;
    
    if noEvent(port_1.m_flow > 0 and port_3.m_flow < 0) then
      port_2.h_a = port_1.h_b;
      port_2.Xi_a = port_1.Xi_b;
    elseif noEvent(port_1.m_flow < 0 and port_3.m_flow > 0) then
      port_2.h_a = port_3.h_b;
      port_2.Xi_a = port_3.Xi_b;
    elseif noEvent(port_1.m_flow > 0 and port_3.m_flow > 0) then
      port_2.h_a = (port_1.m_flow*port_1.h_b + port_3.m_flow*port_3.h_b)/(port_1.m_flow
         + port_3.m_flow);
      port_2.Xi_a = (port_1.m_flow*port_1.Xi_b + port_3.m_flow*port_3.Xi_b)/(port_1.m_flow
         + port_3.m_flow);
    else
      port_2.h_a = port_2.h_b 
        "The theoretical inverse flow direction means in this case that all three flows exit the junction";
      port_2.Xi_a = port_2.Xi_b 
        "The theoretical inverse flow direction means in this case that all three flows exit the junction";
    end if;
    
    if noEvent(port_1.m_flow > 0 and port_2.m_flow < 0) then
      port_3.h_a = port_1.h_b;
      port_3.Xi_a = port_1.Xi_b;
    elseif noEvent(port_1.m_flow < 0 and port_2.m_flow > 0) then
      port_3.h_a = port_2.h_b;
      port_3.Xi_a = port_2.Xi_b;
    elseif noEvent(port_1.m_flow > 0 and port_2.m_flow > 0) then
      port_3.h_a = (port_1.m_flow*port_1.h_b + port_2.m_flow*port_2.h_b)/(port_1.m_flow
         + port_2.m_flow);
      port_3.Xi_a = (port_1.m_flow*port_1.Xi_b + port_2.m_flow*port_2.Xi_b)/(port_1.m_flow
         + port_2.m_flow);
    else
      port_3.h_a = port_3.h_b 
        "The theoretical inverse flow direction means in this case that all three flows exit the junction";
      port_3.Xi_a = port_3.Xi_b 
        "The theoretical inverse flow direction means in this case that all three flows exit the junction";
    end if;
    
  end PartialIdealJunction;
  
  redeclare replaceable partial model extends PartialSource_A 
    "Partial source model with a Port_a" 
    
  equation 
    port.p = medium.p;
    port.h_a = medium.h;
    port.Xi_a = medium.Xi;
    
    // Interface to generic implementation
    port.m_flow = port_m_flow;
  end PartialSource_A;
  
  redeclare replaceable partial model extends PartialSource_B 
    "Partial source model with a Port_b" 
    
  equation 
    port.p = medium.p;
    port.h_b = medium.h;
    port.Xi_b = medium.Xi;
    
    // Interface to generic implementation
    port.m_flow = port_m_flow;
  end PartialSource_B;
  
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
    H_flow[1] = noEvent(port_a.m_flow*(if port_a.m_flow > 0 then port_a.h_a else 
      port_a.h_b));
    H_flow[n + 1] = -noEvent(port_b.m_flow*(if port_b.m_flow > 0 then port_b.h_a else 
            port_b.h_b));
    mXi_flow[1, :] = noEvent(port_a.m_flow*(if port_a.m_flow > 0 then port_a.Xi_a[
      :] else port_a.Xi_b[:]));
    mXi_flow[n + 1, :] = -noEvent(port_b.m_flow*(if port_b.m_flow > 0 then port_b.Xi_a[
      :] else port_b.Xi_b[:]));
    // Quantities related to connectors
    eta_a = if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then 
            eta_nominal else (if use_approxPortProperties then eta[1] else (if 
      m_flow[1] >= 0 then Medium.dynamicViscosity(Medium.setState_phX(
      port_a.p,
      port_a.h_a,
      port_a.Xi_a)) else eta[1])));
    d_a = if use_d_nominal then d_nominal else (if use_approxPortProperties then 
      d[1] else (if m_flow[1] >= 0 then Medium.density_phX(
      port_a.p,
      port_a.h_a,
      port_a.Xi_a) else d[1]));
    eta_b = if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then 
            eta_nominal else (if use_approxPortProperties then eta[n] else (if 
      m_flow[n + 1] < 0 then Medium.dynamicViscosity(Medium.setState_phX(
      port_b.p,
      port_b.h_a,
      port_b.Xi_a)) else eta[n])));
    d_b = if use_d_nominal then d_nominal else (if use_approxPortProperties then 
      d[n] else (if m_flow[n + 1] >= 0 then d[n] else Medium.density_phX(
      port_b.p,
      port_b.h_a,
      port_b.Xi_a)));
    
    // Connector equations
    // Specific enthalpy
    port_a.h_b = medium[1].h;
    port_b.h_b = medium[n].h;
    // Substance mass fractions
    port_a.Xi_b = medium[1].Xi;
    port_b.Xi_b = medium[n].Xi;
    
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
 
</html>",   revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  end PartialSymmetricDistributedPipe;
  
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
    H_flow[1] = noEvent(port_a.m_flow*(if port_a.m_flow > 0 then port_a.h_b else 
      port_a.h_a));
    H_flow[n + 1] = -noEvent(port_b.m_flow*(if port_b.m_flow > 0 then port_b.h_a else 
            port_b.h_b));
    mXi_flow[1, :] = noEvent(port_a.m_flow*(if port_a.m_flow > 0 then port_a.Xi_b[
      :] else port_a.Xi_a[:]));
    mXi_flow[n + 1, :] = -noEvent(port_b.m_flow*(if port_b.m_flow > 0 then port_b.Xi_a[
      :] else port_b.Xi_b[:]));
    // Quantities related to connectors
    eta_a = if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then 
            eta_nominal else (if use_approxPortProperties then eta[1] else (if 
      m_flow[1] >= 0 then Medium.dynamicViscosity(Medium.setState_phX(
      port_a.p,
      port_a.h_b,
      port_a.Xi_b)) else eta[1])));
    d_a = if use_d_nominal then d_nominal else (if use_approxPortProperties then 
      d[1] else (if m_flow[1] >= 0 then Medium.density_phX(
      port_a.p,
      port_a.h_b,
      port_a.Xi_b) else d[1]));
    eta_b = if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then 
            eta_nominal else (if use_approxPortProperties then eta[n] else (if 
      m_flow[n + 1] < 0 then Medium.dynamicViscosity(Medium.setState_phX(
      port_b.p,
      port_b.h_a,
      port_b.Xi_a)) else eta[n])));
    d_b = if use_d_nominal then d_nominal else (if use_approxPortProperties then 
      d[n] else (if m_flow[n + 1] >= 0 then d[n] else Medium.density_phX(
      port_b.p,
      port_b.h_a,
      port_b.Xi_a)));
    
    // Connector equations
    // Specific enthalpy
    port_a.h_a = medium[1].h;
    port_b.h_b = medium[n].h;
    // Substance mass fractions
    port_a.Xi_a = medium[1].Xi;
    port_b.Xi_b = medium[n].Xi;
    
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
 
</html>",   revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  end PartialAsymmetricDistributedPipe;
  
end EffortPairs_A;
