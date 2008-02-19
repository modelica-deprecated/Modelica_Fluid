package TripleEffort_A 
  "Inflow Outflow Outside implementation A (no explicit junctions, asymmetric models, static momentum balances)" 
  extends Interfaces.PartialFluidInterface(usesNewConnectionSemantics=true);
  redeclare replaceable connector extends FluidPort 
    "Interface for quasi one-dimensional fluid flow in a piping network (incompressible or compressible, one or more phases, one or more substances)" 
    
    Medium.AbsolutePressure p "Pressure in the connection point";
    flow Medium.MassFlowRate m_flow 
      "Mass flow rate from the connection point into the component";
    
    // Enthalpy
    Medium.SpecificEnthalpy h_outside 
      "Specific enthalpy outside of component close to port";
  /*   To infuse the symbolically manipulated equations, this was included in the connection semantics
     = noEvent(if m_flow >= 0 then h_outflow else h_inflow) 
*/
    
    /*inside*/
    Medium.SpecificEnthalpy h_inflow 
      "h inside of component if m_flow >= 0, otherwise h_outside";
    
    /*inside*/
    Medium.SpecificEnthalpy h_outflow 
      "h inside of component if m_flow < 0, otherwise h_outside";
    
    flow Medium.EnthalpyFlowRate H_flow 
      "Enthalpy flow rate (positive from port to component)";
  /*  To infuse the symbolically manipulated equations, this was included in the connection semantics
    = semiLinear(m_flow, h_inflow, h_outflow) 
*/
    
    // Composition  
    Medium.MassFraction Xi_outside[Medium.nXi] 
      "Independent mixture mass fractions m_i/m outside of component close to port";
  /*   To infuse the symbolically manipulated equations, this was included in the connection semantics
     = noEvent(if m_flow >= 0 then Xi_outflow else Xi_inflow) 
*/
    
    /*inside*/
    Medium.MassFraction Xi_inflow[Medium.nXi] 
      "Xi inside of component if m_flow >= 0, otherwise Xi_outside";
    
    /*inside*/
    Medium.MassFraction Xi_outflow[Medium.nXi] 
      "Xi inside of component if m_flow < 0, otherwise Xi_outside";
    
    flow Medium.MassFlowRate mXi_flow[Medium.nXi] 
      "Mass flow rates of the independent substances (positive from port to component)";
  /*   To infuse the symbolically manipulated equations, this was included in the connection semantics
     = semiLinear(m_flow, Xi_inflow, Xi_outflow) 
*/
    
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
    "Implements unsupported connection semantics" 
    
  equation 
  /*
  Original form of equations generated from the unsupported connection semantics
  
  // Effort variables _without_ inside prefix: Equality 
  port_a.p = port_b.p;
  port_a.h_outside = port_b.h_outside;
  port_a.Xi_outside = port_b.Xi_outside;
  
  // Effort variables _with_ inside prefix: Ignore
  
  // Flow variables: Sum to zero 
  port_a.m_flow + port_b.m_flow = 0;
  port_a.H_flow + port_b.H_flow = 0;
  port_a.mXi_flow + port_b.mXi_flow = 0;
 
  // Plus declaration equations inside connector (not listed here)
*/
    
    // Manipulated form of connection equations using new connection semantics
    // Taken from "New Modelica_Fluid Principles" by Elmqvist, Casella, Otter, Mattson
    // Has been adapted to the inverted signs of flow variables
    port_a.p = port_b.p;
    port_a.m_flow + port_b.m_flow = 0;
    
    // Enthalpy flow  
    port_a.h_inflow = port_b.h_outflow;
    port_a.h_outflow = port_b.h_inflow;
    
    port_a.h_outside = noEvent(if -port_a.m_flow > 0 then port_a.h_outflow else 
            port_a.h_inflow);
    port_b.h_outside = noEvent(if -port_b.m_flow > 0 then port_b.h_outflow else 
            port_b.h_inflow);
    
    port_a.H_flow = -semiLinear(
            -port_a.m_flow,
            port_a.h_inflow,
            port_a.h_outflow);
    port_b.H_flow = -semiLinear(
            -port_b.m_flow,
            port_b.h_inflow,
            port_b.h_outflow);
    
    // Substance mass flow
    port_a.Xi_inflow = port_b.Xi_outflow;
    port_a.Xi_outflow = port_b.Xi_inflow;
    
    port_a.Xi_outside = noEvent(if -port_a.m_flow > 0 then port_a.Xi_outflow else 
            port_a.Xi_inflow);
    port_b.Xi_outside = noEvent(if -port_b.m_flow > 0 then port_b.Xi_outflow else 
            port_b.Xi_inflow);
    
    port_a.mXi_flow = -semiLinear(
            -port_a.m_flow,
            port_a.Xi_inflow,
            port_a.Xi_outflow);
    port_b.mXi_flow = -semiLinear(
            -port_b.m_flow,
            port_b.Xi_inflow,
            port_b.Xi_outflow);
    
    annotation (
      defaultComponentName="semantics",
      Coordsys(extent=[-100,-20; 100,20], scale=0.05),
      Icon(Line(points=[-90,0; 90,0], style(
            color=69,
            rgbcolor={0,128,255},
            thickness=2)), Text(extent=[-300,0; 300,-64], string=
              "new semantics")),
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
    port_a.h_outflow = medium.h*ones(n_a);
    port_b.h_outflow = medium.h*ones(n_b);
    
    // Substance mass flows
    for i in 1:n_a loop
      port_a[i].Xi_outflow = medium.Xi;
    end for;
    for i in 1:n_b loop
      port_b[i].Xi_outflow = medium.Xi;
    end for;
    
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
    // Mass balance
    port_a.m_flow + port_b.m_flow = 0;
    
    // Energy balance when fluid flows from port_a to port_b
    port_b.h_outflow = port_a.h_inflow;
    // Energy balance when fluid flows from port_b to port_a
    port_a.h_outflow = port_b.h_inflow;
    
    // Substance mass balance when fluid flows from port_a to port_b
    port_b.Xi_outflow = port_a.Xi_inflow;
    // Substance mass balance when fluid flows from port_b to port_a
    port_a.Xi_outflow = port_b.Xi_inflow;
    
    // Design direction of mass flow rate
    m_flow = port_a.m_flow;
    
    // Pressure difference between ports
    dp = port_a.p - port_b.p;
    
    // This approach provides upstream and downstream properties
    p_designDirection = port_a.p 
      "Upstream pressure if flow is in design direction";
    h_designDirection = port_a.h_inflow 
      "Upstream specific enthalpy if flow is in design direction";
    Xi_designDirection = port_a.Xi_inflow 
      "Upstream mass fractions if flow is in design direction";
    p_nonDesignDirection = port_b.p 
      "Upstream pressure if flow is in non-design direction";
    h_nonDesignDirection = port_b.h_inflow 
      "Upstream specific enthalpy if flow is in non-design direction";
    Xi_nonDesignDirection = port_b.Xi_inflow 
      "Upstream mass fractions if flow is in non-design direction";
    
  end PartialTransportIsenthalpic;
  
  redeclare replaceable partial model extends PartialTransportIsenthalpicAA 
    "Partial isenthalpic element transporting fluid between two ports without storing mass or energy (two Port_a's, not yet implemented)" 
    
  equation 
    assert(false, "The PartialTransportIsenthalpicAA was not yet implemented for this approach.");
    
    annotation (Icon(Rectangle(extent=[-102,102; 102,-102], style(
            color=1,
            rgbcolor={255,0,0},
            pattern=2,
            thickness=2))));
  end PartialTransportIsenthalpicAA;
  
  redeclare replaceable partial model extends PartialTransportIsenthalpicAB 
    "Partial isenthalpic element transporting fluid between two ports without storing mass or energy (a Port_a and Port_b each, not yet implemented)" 
    
  equation 
    assert(false, "The PartialTransportIsenthalpicAB was not yet implemented for this approach.");
    
    annotation (Icon(Rectangle(extent=[-102,102; 102,-102], style(
            color=1,
            rgbcolor={255,0,0},
            pattern=2,
            thickness=2))));
  end PartialTransportIsenthalpicAB;
  
  redeclare replaceable partial model extends PartialIdealJunction 
    "Partial infinitesimal junction model" 
    
  equation 
    assert(false, "The PartialIdealJunction was not yet implemented for this approach.");
    
    annotation (Icon(Rectangle(extent=[-102,102; 102,-102], style(
            color=1,
            rgbcolor={255,0,0},
            pattern=2,
            thickness=2))));
  end PartialIdealJunction;
  
  redeclare replaceable partial model extends PartialSource_A 
    "Partial source model with a Port_a" 
    
  equation 
    port.p = medium.p;
    port.h_outflow = medium.h;
    port.Xi_outflow = medium.Xi;
    
    // Interface to generic implementation
    port.m_flow = port_m_flow;
  end PartialSource_A;
  
  redeclare replaceable partial model extends PartialSource_B 
    "Partial source model with a Port_b" 
    
  equation 
    port.p = medium.p;
    port.h_outflow = medium.h;
    port.Xi_outflow = medium.Xi;
    
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
    eta_a = if not WallFriction.use_eta then 1.e-10 else (if 
      use_eta_nominal then eta_nominal else (if use_approxPortProperties then 
            eta[1] else (if m_flow[1] >= 0 then Medium.dynamicViscosity(
      Medium.setState_phX(
            port_a.p,
            port_a.h_inflow,
            port_a.Xi_inflow)) else eta[1])));
    eta_b = if not WallFriction.use_eta then 1.e-10 else (if 
      use_eta_nominal then eta_nominal else (if use_approxPortProperties then 
            eta[n] else (if m_flow[n + 1] < 0 then Medium.dynamicViscosity(
      Medium.setState_phX(
            port_b.p,
            port_b.h_inflow,
            port_b.Xi_inflow)) else eta[n])));
    d_a = if use_d_nominal then d_nominal else (if use_approxPortProperties then 
            d[1] else (if m_flow[1] >= 0 then Medium.density_phX(
            port_a.p,
            port_a.h_inflow,
            port_a.Xi_inflow) else d[1]));
    d_b = if use_d_nominal then d_nominal else (if use_approxPortProperties then 
            d[n] else (if m_flow[n + 1] >= 0 then d[n] else 
      Medium.density_phX(
            port_b.p,
            port_b.h_inflow,
            port_b.Xi_inflow)));
    
    // Connector equations
    // Enthalpy flow
    port_a.h_outflow = medium[1].h;
    port_b.h_outflow = medium[n].h;
    // Substance mass flow
    port_a.Xi_outflow = medium[1].Xi;
    port_b.Xi_outflow = medium[n].Xi;
    
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
    eta_a = if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then 
            eta_nominal else (if use_approxPortProperties then eta[1] else (if 
      m_flow[1] >= 0 then Medium.dynamicViscosity(Medium.setState_phX(
      port_a.p,
      port_a.h_inflow,
      port_a.Xi_inflow)) else eta[1])));
    eta_b = if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then 
            eta_nominal else (if use_approxPortProperties then eta[n] else (if 
      m_flow[n + 1] < 0 then Medium.dynamicViscosity(Medium.setState_phX(
      port_b.p,
      port_b.h_inflow,
      port_b.Xi_inflow)) else eta[n])));
    d_a = if use_d_nominal then d_nominal else (if use_approxPortProperties then 
      d[1] else (if m_flow[1] >= 0 then Medium.density_phX(
      port_a.p,
      port_a.h_inflow,
      port_a.Xi_inflow) else d[1]));
    d_b = if use_d_nominal then d_nominal else (if use_approxPortProperties then 
      d[n] else (if m_flow[n + 1] >= 0 then d[n] else Medium.density_phX(
      port_b.p,
      port_b.h_inflow,
      port_b.Xi_inflow)));
    
  // Connector equations
  // Enthalpy flow
    port_a.h_outflow = medium[1].h;
    port_b.h_outflow = medium[n].h;
  // Substance mass flow
    port_a.Xi_outflow = medium[1].Xi;
    port_b.Xi_outflow = medium[n].Xi;
    
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
end TripleEffort_A;
