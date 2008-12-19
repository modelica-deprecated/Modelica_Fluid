within Modelica_Fluid;
package Pipes "Devices for conveying fluid"
    extends Modelica_Fluid.Icons.VariantLibrary;

  model StaticPipe
    "Basic pipe flow model without storage of momentum, mass or energy"

    // extending PartialStraightPipe
    extends Modelica_Fluid.Pipes.BaseClasses.PartialStraightPipe(
          redeclare model HeatTransfer = 
          BaseClasses.HeatTransfer.NoHeatTransfer (                            nParallel=nParallel));

    // Assumptions
    parameter Types.Dynamics momentumDynamics=
      system.momentumDynamics
      "Formulation of momentum balance, if pressureLoss options available" 
      annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));

    // Initialization
    parameter Medium.AbsolutePressure p_a_start=system.p_start
      "Start value of pressure at port a" 
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_b_start=p_a_start
      "Start value of pressure at port b" 
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.MassFlowRate m_flow_start = system.m_flow_start
      "Start value for mass flow rate" 
       annotation(Evaluate=true, Dialog(tab = "Initialization"));

    PressureLoss pressureLoss(
            redeclare final package Medium = Medium,
            final n=1,
            state={Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)),
                   Medium.setState_phX(port_b.p, inStream(port_b.h_outflow), inStream(port_b.Xi_outflow))},
            final allowFlowReversal=allowFlowReversal,
            final momentumDynamics=momentumDynamics,
            final p_a_start=p_a_start,
            final p_b_start=p_b_start,
            final m_flow_start=m_flow_start,
            final nParallel=nParallel,
            final length={length},
            final crossArea={crossArea, crossArea},
            final perimeter={perimeter},
            final roughness={roughness},
            final height_ab=height_ab,
            final g=system.g) "Pressure loss model" 
       annotation (Placement(transformation(extent={{-38,-18},{38,18}},rotation=0)));
  equation
    port_a.m_flow = pressureLoss.m_flow[1];
    0 = port_a.m_flow + port_b.m_flow;
    port_a.h_outflow = inStream(port_b.h_outflow);
    port_b.h_outflow = inStream(port_a.h_outflow);
    port_a.Xi_outflow = inStream(port_b.Xi_outflow);
    port_b.Xi_outflow = inStream(port_a.Xi_outflow);
    port_a.C_outflow = inStream(port_b.C_outflow);
    port_b.C_outflow = inStream(port_a.C_outflow);
    annotation (defaultComponentName="pipe", Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
              -100},{100,100}}),      graphics));
  end StaticPipe;

  model LumpedPipe "Example for a composite pipe model"

    // extending PartialStraightPipe
    extends Modelica_Fluid.Pipes.BaseClasses.PartialStraightPipe;

    // Assumptions
    parameter Types.Dynamics energyDynamics=system.energyDynamics
      "Formulation of energy balance" 
      annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
    parameter Types.Dynamics massDynamics=system.massDynamics
      "Formulation of mass balance" 
      annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
    parameter Types.Dynamics momentumDynamics=system.momentumDynamics
      "Formulation of momentum balance, if pressureLoss options available" 
      annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));

    // Initialization
    parameter Medium.AbsolutePressure p_a_start=system.p_start
      "Start value of pressure at port a" 
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_b_start=p_a_start
      "Start value of pressure at port b" 
      annotation(Dialog(tab = "Initialization"));

    parameter Boolean use_T_start=true "Use T_start if true, otherwise h_start"
       annotation(Evaluate=true, Dialog(tab = "Initialization"));
    parameter Medium.Temperature T_start=if use_T_start then system.T_start else 
                Medium.temperature_phX(
          (p_a_start + p_b_start)/2,
          h_start,
          X_start) "Start value of temperature" 
      annotation(Evaluate=true, Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=if use_T_start then 
          Medium.specificEnthalpy_pTX(
          (p_a_start + p_b_start)/2,
          T_start,
          X_start) else Medium.h_default "Start value of specific enthalpy" 
      annotation(Evaluate=true, Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX]=Medium.X_default
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    parameter Medium.ExtraProperty C_start[Medium.nC](
         quantity=Medium.extraPropertiesNames)=fill(0, Medium.nC)
      "Start value of trace substances" 
      annotation (Dialog(tab="Initialization", enable=Medium.nC > 0));

    parameter Medium.MassFlowRate m_flow_start = system.m_flow_start
      "Start value for mass flow rate" 
       annotation(Evaluate=true, Dialog(tab = "Initialization"));

    HeatTransfer heatTransfer(
      redeclare final package Medium = Medium,
      final n=1,
      final nParallel=nParallel,
      final length={length},
      final crossArea={crossArea, crossArea},
      final perimeter={perimeter},
      final roughness={roughness},
      final use_fluidHeatPort=true,
      state={volume.medium.state},
      m_flow = {0.5*(port_a.m_flow - port_b.m_flow)}) "Heat transfer model" 
        annotation (Placement(transformation(extent={{-11,14},{11,36}},rotation=0)));

    StaticPipe staticPipe1(
      redeclare package Medium = Medium,
      allowFlowReversal=allowFlowReversal,
      momentumDynamics=momentumDynamics,
      nParallel=nParallel,
      length=length/2,
      roughness=roughness,
      diameter=diameter,
      perimeter=perimeter,
      crossArea=crossArea,
      height_ab=height_ab/2,
      m_flow_start=m_flow_start,
      redeclare final model PressureLoss = PressureLoss) 
      annotation (Placement(transformation(extent={{-60,-40},{-40,-20}},
            rotation=0)));
    Modelica_Fluid.Vessels.Volume volume(
      redeclare package Medium = Medium,
      energyDynamics=energyDynamics,
      massDynamics=massDynamics,
      p_start=(p_a_start+p_b_start)/2,
      use_T_start=use_T_start,
      T_start=T_start,
      h_start=h_start,
      X_start=X_start,
      C_start=C_start,
      V=V,
      nPorts=2,
      portDiameters={0,0},
      neglectPortDiameters=true) 
      annotation (Placement(transformation(extent={{-10,-20},{10,0}},  rotation=
             0)));
    StaticPipe staticPipe2(
      redeclare package Medium = Medium,
      allowFlowReversal=allowFlowReversal,
      momentumDynamics=momentumDynamics,
      nParallel=nParallel,
      length=length/2,
      roughness=roughness,
      diameter=diameter,
      perimeter=perimeter,
      crossArea=crossArea,
      height_ab=height_ab/2,
      m_flow_start=m_flow_start,
      redeclare final model PressureLoss = PressureLoss)   annotation (Placement(transformation(extent={{40,-40},
              {60,-20}},         rotation=0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort 
      annotation (Placement(transformation(extent={{-10,44},{10,64}}, rotation=
              0)));
  equation
    connect(staticPipe1.port_a, port_a) 
      annotation (Line(points={{-60,-30},{-80,-30},{-80,0},{-100,0}},
                                                  color={0,127,255}));
    connect(staticPipe2.port_b, port_b) 
      annotation (Line(points={{60,-30},{80,-30},{80,0},{100,0}},
                                                color={0,127,255}));
    connect(heatPort, heatTransfer.wallHeatPort[1]) annotation (Line(
        points={{0,54},{0,32.7}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(heatTransfer.fluidHeatPort[1], volume.heatPort) annotation (Line(
        points={{0,18.4},{0,-0.2}},
        color={191,0,0},
        smooth=Smooth.None));
    annotation (defaultComponentName="pipe",Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={Ellipse(
            extent={{-10,10},{10,-10}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid)}),Documentation(info="<html>
<p>
Simple pipe model consisting of one volume, 
wall friction (with different friction correlations)
and gravity effect. This model is mostly used to demonstrate how
to build up more detailed models from the basic components.
Note, if the \"heatPort\" is not connected, then the pipe
is totally insulated (= no thermal flow from the fluid to the
pipe wall/environment).
</p>
</html>"),
      Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics));
    connect(staticPipe1.port_b, volume.ports[1])   annotation (Line(
        points={{-40,-30},{0,-30},{0,-18}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(staticPipe2.port_a, volume.ports[2])   annotation (Line(
        points={{40,-30},{0,-30},{0,-22}},
        color={0,127,255},
        smooth=Smooth.None));
  end LumpedPipe;

  model DistributedPipe "Distributed pipe model"

    import Modelica_Fluid.Types;
    import Modelica_Fluid.Types.ModelStructure;

    // extending PartialStraightPipe
    extends Modelica_Fluid.Pipes.BaseClasses.PartialStraightPipe(
      final port_a_exposesState = (modelStructure == ModelStructure.av_b) or (modelStructure == ModelStructure.av_vb),
      final port_b_exposesState = (modelStructure == ModelStructure.a_vb) or (modelStructure == ModelStructure.av_vb));

    // distributed volume model
    extends Modelica_Fluid.Vessels.BaseClasses.PartialDistributedVolume(
      final n = nNodes,
      Qs_flow = heatTransfer.Q_flow);

    // Assumptions
    parameter Types.Dynamics momentumDynamics=system.momentumDynamics
      "Formulation of momentum balances, if pressureLoss options available" 
      annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));

    // Initialization
    parameter Medium.MassFlowRate m_flow_start = system.m_flow_start
      "Start value for mass flow rate" 
       annotation(Evaluate=true, Dialog(tab = "Initialization"));

    // Discretization
    parameter Integer nNodes(min=1)=2 "Number of discrete flow volumes" 
      annotation(Dialog(tab="Advanced"),Evaluate=true);

    parameter Types.ModelStructure modelStructure=Types.ModelStructure.av_vb
      "Determines whether flow or volume models are present at the ports" 
      annotation(Dialog(tab="Advanced"), Evaluate=true);

    parameter Boolean lumpedPressure=false
      "=true to lump all pressure states into one" 
      annotation(Dialog(tab="Advanced"),Evaluate=true);
    final parameter Integer nFlows=if lumpedPressure then nFlowsLumped else nFlowsDistributed
      "number of flow models in pressureLoss";
    final parameter Integer nFlowsDistributed=if modelStructure==Types.ModelStructure.a_v_b then n+1 else if (modelStructure==Types.ModelStructure.a_vb or modelStructure==Types.ModelStructure.av_b) then n else n-1;
    final parameter Integer nFlowsLumped=if modelStructure==Types.ModelStructure.a_v_b then 2 else 1;
    final parameter Integer iLumped=integer(n/2)+1
      "Index of control volume with representative state if lumpedPressure" 
      annotation(Evaluate=true);

    // Advanced model options
    parameter Boolean useInnerPortProperties=false
      "=true to take port properties for pressure drops from internal control volumes"
      annotation(Dialog(tab="Advanced"),Evaluate=true);
    Medium.ThermodynamicState state_a "state defined by volume outside port_a";
    Medium.ThermodynamicState state_b "state defined by volume outside port_b";
    Medium.ThermodynamicState[nFlows+1] flowState
      "state vector for pressureLoss model";

    // Pressure loss model
    PressureLoss pressureLoss(
            redeclare final package Medium = Medium,
            final n=nFlows,
            state=flowState,
            final allowFlowReversal=allowFlowReversal,
            final momentumDynamics=momentumDynamics,
            final p_a_start=p_a_start,
            final p_b_start=p_b_start,
            final m_flow_start=m_flow_start,
            final nParallel=nParallel,
            final length=if modelStructure == Types.ModelStructure.a_v_b then 
                              cat(1, {length/n/2}, fill(length/n, n-1), {length/n/2}) else 
                              fill(length/nFlows, nFlows),
            final crossArea=fill(crossArea, nFlows+1),
            final perimeter=fill(perimeter, nFlows),
            final roughness=fill(roughness, nFlows),
            final height_ab=height_ab,
            final g=system.g) "Pressure loss model" 
       annotation (Placement(transformation(extent={{-77,-57},{77,-23}},rotation=0)));

    // Flow quantities
    Medium.MassFlowRate[n+1] m_flow(
       each min=if allowFlowReversal then -Modelica.Constants.inf else 0,
       each start=m_flow_start)
      "Mass flow rates of fluid across segment boundaries";
    Medium.MassFlowRate[n+1, Medium.nXi] mXi_flow
      "Independent mass flow rates across segment boundaries";
    Medium.MassFlowRate[n+1, Medium.nC] mC_flow
      "Trace substance mass flow rates across segment boundaries";
    Medium.EnthalpyFlowRate[n+1] H_flow
      "Enthalpy flow rates of fluid across segment boundaries";

    // Wall heat transfer
    Interfaces.HeatPorts_a[nNodes] heatPorts 
      annotation (Placement(transformation(extent={{-10,44},{10,64}}), iconTransformation(extent={{-30,44},{32,60}})));

    HeatTransfer heatTransfer(
      redeclare each final package Medium = Medium,
      final n=n,
      final nParallel=nParallel,
      final length=fill(length/n, n),
      final crossArea=fill(crossArea, n+1),
      final perimeter=fill(perimeter, n),
      final roughness=fill(roughness, n),
      state=medium.state,
      m_flow = 0.5*(m_flow[1:n]+m_flow[2:n+1])) "Heat transfer model" 
        annotation (Placement(transformation(extent={{-20,-5},{20,35}},  rotation=0)));

  equation
    assert(nNodes > 1 or modelStructure <> ModelStructure.av_vb,
       "nNodes needs to be at least 2 for modelStructure av_vb, as flow model disappears otherwise!");

    // distributed volume
    if modelStructure == Types.ModelStructure.av_vb then
      // half volumes at the ports for modelStructure av_vb
      fluidVolume = cat(1, {V/(n-1)/2}, fill(V/(n-1), n-2), {V/(n-1)/2});
    else
      // even distribution of volumes else
      fluidVolume=fill(V/n, n);
    end if;
    // Source/sink terms for mass and energy balances
    Ws_flow=zeros(n);
    for i in 1:n loop
      ms_flow[i] = m_flow[i] - m_flow[i + 1];
      msXi_flow[i, :] = mXi_flow[i, :] - mXi_flow[i + 1, :];
      msC_flow[i, :]  = mC_flow[i, :]  - mC_flow[i + 1, :];
      Hs_flow[i] = H_flow[i] - H_flow[i + 1];
    end for;

    // Distributed flow quantities, upwind discretization
    for i in 2:n loop
      H_flow[i] = semiLinear(m_flow[i], medium[i - 1].h, medium[i].h);
      mXi_flow[i, :] = semiLinear(m_flow[i], medium[i - 1].Xi, medium[i].Xi);
      mC_flow[i, :]  = semiLinear(m_flow[i], C[i - 1, :],         C[i, :]);
    end for;
    H_flow[1] = semiLinear(port_a.m_flow, inStream(port_a.h_outflow), medium[1].h);
    H_flow[n + 1] = -semiLinear(port_b.m_flow, inStream(port_b.h_outflow), medium[n].h);
    mXi_flow[1, :] = semiLinear(port_a.m_flow, inStream(port_a.Xi_outflow), medium[1].Xi);
    mXi_flow[n + 1, :] = -semiLinear(port_b.m_flow, inStream(port_b.Xi_outflow), medium[n].Xi);
    mC_flow[1, :] = semiLinear(port_a.m_flow, inStream(port_a.C_outflow), C[1, :]);
    mC_flow[n + 1, :] = -semiLinear(port_b.m_flow, inStream(port_b.C_outflow), C[n, :]);

    // Boundary conditions
    port_a.m_flow    = m_flow[1];
    port_b.m_flow    = -m_flow[n + 1];
    port_a.h_outflow = medium[1].h;
    port_b.h_outflow = medium[n].h;
    port_a.Xi_outflow = medium[1].Xi;
    port_b.Xi_outflow = medium[n].Xi;
    port_a.C_outflow = C[1, :];
    port_b.C_outflow = C[n, :];
    // The two equations below are not correct if C is stored in volumes.
    // C should be treated the same way as Xi.
    //port_a.C_outflow = inStream(port_b.C_outflow);
    //port_b.C_outflow = inStream(port_a.C_outflow);

    if useInnerPortProperties and n > 0 then
      state_a = Medium.setState_phX(port_a.p, medium[1].h, medium[1].Xi);
      state_b = Medium.setState_phX(port_b.p, medium[n].h, medium[n].Xi);
    else
      state_a = Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow));
      state_b = Medium.setState_phX(port_b.p, inStream(port_b.h_outflow), inStream(port_b.Xi_outflow));
    end if;

    if lumpedPressure then
      if modelStructure <> ModelStructure.av_vb then
        // all pressures are equal
        fill(medium[1].p, n-1) = medium[2:n].p;
      elseif n > 2 then
        // need two pressures
        fill(medium[1].p, iLumped-2) = medium[2:iLumped-1].p;
        fill(medium[n].p, n-iLumped) = medium[iLumped:n-1].p;
      end if;
      if modelStructure == ModelStructure.a_v_b then
        m_flow[1] = pressureLoss.m_flow[1];
        flowState[1] = state_a;
        flowState[2] = medium[iLumped].state;
        flowState[3] = state_b;
        m_flow[n+1] = pressureLoss.m_flow[2];
      elseif modelStructure == ModelStructure.av_b then
        port_a.p = medium[1].p;
        flowState[1] = medium[iLumped].state;
        flowState[2] = state_b;
        m_flow[n+1] = pressureLoss.m_flow[1];
      elseif modelStructure == ModelStructure.a_vb then
        m_flow[1] = pressureLoss.m_flow[1];
        flowState[1] = state_a;
        flowState[2] = medium[iLumped].state;
        port_b.p = medium[n].p;
      else // av_vb
        port_a.p = medium[1].p;
        flowState[1] = medium[1].state;
        m_flow[iLumped] = pressureLoss.m_flow[1];
        flowState[2] = medium[n].state;
        port_b.p = medium[n].p;
      end if;
    else
      if modelStructure == ModelStructure.a_v_b then
        flowState[1] = state_a;
        flowState[2:n+1] = medium[1:n].state;
        flowState[n+2] = state_b;
        //m_flow = pressureLoss.m_flow;
        for i in 1:n+1 loop
          m_flow[i] = pressureLoss.m_flow[i];
        end for;
      elseif modelStructure == ModelStructure.av_b then
        flowState[1:n] = medium[1:n].state;
        flowState[n+1] = state_b;
        //m_flow[2:n+1] = pressureLoss.m_flow;
        for i in 2:n+1 loop
          m_flow[i] = pressureLoss.m_flow[i-1];
        end for;
        port_a.p = medium[1].p;
      elseif modelStructure == ModelStructure.a_vb then
        flowState[1] = state_a;
        flowState[2:n+1] = medium[1:n].state;
        //m_flow[1:n] = pressureLoss.m_flow;
        for i in 1:n loop
          m_flow[i] = pressureLoss.m_flow[i];
        end for;
        port_b.p = medium[n].p;
      else // av_vb
        flowState[1:n] = medium[1:n].state;
        //m_flow[2:n] = pressureLoss.m_flow[1:n-1];
        for i in 2:n loop
          m_flow[i] = pressureLoss.m_flow[i-1];
        end for;
        port_a.p = medium[1].p;
        port_b.p = medium[n].p;
      end if;
    end if;

    connect(heatPorts, heatTransfer.wallHeatPort) 
      annotation (Line(points={{0,54},{0,29}}, color={191,0,0}));
    annotation (defaultComponentName="pipe",
  Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={Ellipse(
            extent={{-72,10},{-52,-10}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid), Ellipse(
            extent={{50,10},{70,-10}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid)}),
  Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{100,
              100}},
          grid={1,1}),
          graphics),
  Documentation(info="<html>
<p>Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PartialStraightPipe\">PartialStraightPipe</a>. 
The total volume is determined by geometry parameters. It is split into nNodes pipe segments of equal size along the flow path. 
The default value is nNodes=2.
<p><b>Mass and Energy balance</b></p>
One mass and one energy balance if formulated for each pipe segment. 
The mass and energy balances are inherited from <a href=\"Modelica:Modelica_Fluid.Vessels.BaseClasses.PartialDistributedVolume\">PartialDistributedVolume</a>. 
The additional component <b><tt>HeatTransfer</tt></b> specifies the source term <tt>Qs_flow</tt> in the energy balance. 
The default component uses a constant coefficient for the heat transfer between the bulk flow and the segment boundaries exposed through the <tt>heatPorts</tt>. 
The <tt>HeatTransfer</tt> model is replaceable and can be exchanged with any model extended from 
<a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialHeatTransfer\">BaseClasses.PartialHeatTransfer</a>.</p>
<p><b>Momentum balance</b></p>
The momentum balance is determined by the <b><tt>PressureLoss</tt></b> component, which can be replaced with any model extended from 
<a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PressureLoss.PartialPressureLoss\">BaseClasses.PartialPressureLoss</a>.
The default setting is steady-state <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PressureLoss.DetailedFlow\">DetailedFlow</a>.
The momentum balances are formed across the segment boundaries along the flow path according to the staggered grid approach. 
The default symmetric model is characterized by one momentum balance inside the pipe with nNodes=2 fluid segments.
An alternative symmetric variation with nNodes+1 momentum balances, half a momentum balance at each port, as well as  
non-symmetric variations can be obtained by chosing a different value for the parameter <tt><b>modelStructure</b></tt>. 
The options include:
<ul>
<li><tt>av_vb</tt>: nNodes-1 equal momentum balances between nNodes pipe segments. 
This results in potential pressure states at both ports.
<li><tt>a_v_b</tt>: Alternative symmetric setting with nNodes+1 momentum balances across nNodes equal pipe segments, 
half a momentum balance at each port. Connecting two pipes therefore results in algebraic pressures at the ports. 
The specification of good start values for the port pressures is essential in order to solve large systems.</li>
<li><tt>av_b</tt>: nNodes momentum balances, one between nth volume and <tt>port_b</tt>, potential pressure state at <tt>port_a</tt></li>
<li><tt>a_vb</tt>: nNodes momentum balance, one between first volume and <tt>port_a</tt>, potential pressure state at <tt>port_b</tt></li>
</ul></p>
 
<p>The PressureLoss contains
<ul>
<li>pressure drop due to friction and other dissipative losses</li>
<li>gravity effects for non-horizontal pipes</li>
</ul>
It does not model changes in pressure resulting from significant variation of flow velocity along the flow path (with the assumption of a constant cross sectional area it must result from fluid density changes, such as in two-phase flow).
 
When connecting two components, e.g. two pipes, the momentum balance across the connection point reduces to
</p> 
<pre>pipe1.port_b.p = pipe2.port_a.p</pre>
<p>
This is only true if the flow velocity remains the same on each side of the connection. 
Consider using a fitting, like <a href=\"Modelica:Modelica_Fluid.Fittings.SuddenExpansion\">SuddenExpansion</a>
for any significant change in diameter, if the resulting effects, such as change in kinetic energy, cannot be neglected. 
This also allows for taking into account friction losses with respect to the actual geometry of the connection point.
</p>
 
</html>",
      revisions="<html>
<ul>
<li><i>5 Dec 2008</i>
    by Michael Wetter:<br>
       Modified mass balance for trace substances. With the new formulation, the trace substances masses <tt>mC</tt> are stored
       in the same way as the species <tt>mXi</tt>.</li>
<li><i>4 Dec 2008</i>
    by R&uuml;diger Franke:<br>
       Derived model from original DistributedPipe models</li>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));

  end DistributedPipe;

  package BaseClasses
    extends Modelica_Fluid.Icons.BaseClassLibrary;

    partial model PartialStraightPipe "Base class for straight pipe models"
      extends Modelica_Fluid.Interfaces.PartialTwoPort;

      // Geometry

      // Note: define nParallel as Real to support inverse calculations
      parameter Real nParallel(min=1)=1 "Number of identical parallel pipes" 
        annotation(Dialog(group="Geometry"));
      parameter SI.Length length "Length" 
        annotation(Dialog(tab="General", group="Geometry"));
      parameter Boolean isCircular=true
        "= true if cross sectional area is circular" 
        annotation (Evaluate, Dialog(tab="General", group="Geometry"));
      parameter SI.Diameter diameter "Diameter of circular pipe" 
        annotation(Dialog(group="Geometry", enable=isCircular));
      parameter SI.Area crossArea=Modelica.Constants.pi*diameter*diameter/4
        "Inner cross section area" 
        annotation(Dialog(tab="General", group="Geometry", enable=not isCircular));
      parameter SI.Length perimeter=Modelica.Constants.pi*diameter
        "Inner perimeter" 
        annotation(Dialog(tab="General", group="Geometry", enable=not isCircular));
      parameter SI.Height roughness(min=0)=2.5e-5
        "Average height of surface asperities (default = smooth steel pipe)" 
          annotation(Dialog(group="Geometry"));
      final parameter SI.Volume V=crossArea*length*nParallel "volume size";

      // Static head
      parameter SI.Length height_ab=0 "Height(port_b) - Height(port_a)" 
          annotation(Dialog(group="Static head"), Evaluate=true);

      // Pressure loss
      replaceable model PressureLoss = 
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.WallFrictionPressureLoss 
        constrainedby BaseClasses.PressureLoss.PartialPressureLoss
        "Characteristics of wall friction and gravity" 
          annotation(Dialog(group="Pressure loss"), choicesAllMatching=true);

      // Heat transfer
      replaceable model HeatTransfer = 
          Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.ConstantHeatTransfer 
        constrainedby
        Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialHeatTransfer
        "Wall heat transfer" 
          annotation (Dialog(group="Heat transfer"),choicesAllMatching=true);
    equation
      assert(length >= height_ab, "Parameter length must be greater or equal height_ab.");

      annotation (defaultComponentName="pipe",Icon(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={1,1}), graphics={
            Rectangle(
              extent={{-100,46},{100,-47}},
              lineColor={0,0,0},
              fillPattern=FillPattern.HorizontalCylinder,
              fillColor={192,192,192}),
            Rectangle(
              extent={{-100,40},{100,-40}},
              lineColor={0,0,0},
              fillPattern=FillPattern.HorizontalCylinder,
              fillColor={0,127,255}),
            Text(
              extent={{-150,-92},{150,-132}},
              lineColor={0,0,255},
              fillPattern=FillPattern.HorizontalCylinder,
              fillColor={0,127,255},
              textString="%name")}),        Documentation(info="<html>
<p>
Base class for one dimensional flow models. It specializes a PartialTwoPort with a parameter interface and icon graphics.
</p>
</html>"),
        Diagram(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={1,1}), graphics));

    end PartialStraightPipe;

    package PressureLoss
      "Pressure loss models for pipes, including wall friction and static head"
          partial model PartialPressureLoss
        "Base class for pipe wall friction models"

            replaceable package Medium = 
            Modelica.Media.Interfaces.PartialMedium "fluid medium" 
              annotation(Dialog(tab="Internal Interface", enable=false));

            input Medium.ThermodynamicState[n+1] state
          "states along design flow";
            output Medium.MassFlowRate[n] m_flow(each start = m_flow_start)
          "mass flow rates along design flow";

            // Discretization
            parameter Integer n=1 "number of flow segments" 
              annotation(Dialog(tab="Internal Interface", enable=false));

            // Mandadory geometry parameters
            parameter Real nParallel "number of parallel pipes" 
               annotation(Dialog(tab="Internal Interface", enable=false,group="Geometry"));
            parameter SI.Length[n] length "Length of segments along flow path" 
               annotation(Dialog(tab="Internal Interface", enable=false,group="Geometry"));
            parameter SI.Area[n+1] crossArea
          "Cross flow area at segment boundaries" 
              annotation(Dialog(tab="Internal Interface", enable=false,group="Geometry"));
            parameter SI.Length[n] perimeter
          "Mean perimeter of segments, used for hydraulic diameter" 
              annotation(Dialog(tab="Internal Interface", enable=false,group="Geometry"));
            parameter SI.Height[n] roughness(each min=0)
          "Average height of surface asperities" 
                annotation(Dialog(tab="Internal Interface", enable=false,group="Geometry",enable=WallFriction.use_roughness));
            parameter SI.Length height_ab
          "Height(state[n+1]) - Height(state[1])" 
                annotation(Dialog(tab="Internal Interface", enable=false,group="Static head"));

            // Additional parameters
            // Note: no outer system is used for default values,
            // as a PressureLoss model is intended as sub-component of other models
            parameter SI.Acceleration g "Constant gravity acceleration" 
              annotation(Dialog(tab="Internal Interface", enable=false,group="Static head"));
            parameter Boolean allowFlowReversal
          "= true to allow flow reversal, false restricts to design direction (state[1] -> state[n+1])"
              annotation(Dialog(tab="Internal Interface", enable=false,group="Assumptions"), Evaluate=true);
            parameter Modelica_Fluid.Types.Dynamics momentumDynamics
          "Formulation of momentum balance, if options available" 
              annotation(Dialog(tab="Internal Interface", enable=false,group = "Assumptions"), Evaluate=true);
            parameter Medium.AbsolutePressure p_a_start
          "Start value for p[1] at design inflow" 
              annotation(Dialog(tab="Internal Interface", enable=false,group = "Initialization"));
            parameter Medium.AbsolutePressure p_b_start
          "Start value for p[n+1] at design outflow" 
              annotation(Dialog(tab="Internal Interface", enable=false,group = "Initialization"));
            parameter Medium.MassFlowRate m_flow_start
          "Start value of mass flow rate" 
              annotation(Dialog(tab="Internal Interface", enable=false,group = "Initialization"));

            // Variables
            //final parameter Medium.AbsolutePressure[n+1] p_start = if n > 0 then linspace(p_a_start, p_b_start, n+1) else {(p_a_start + p_b_start)/2};
            // Note: don't use start values for p to get same behavior as with former PressureLosses.WallFrictionAndGravity
            Medium.AbsolutePressure[n+1] p "pressures of states";
            Medium.AbsolutePressure[n] dp(each start = (p_a_start - p_b_start)/n)
          "pressure loss between states";
          equation
            p = Medium.pressure(state);
            dp = p[1:n] - p[2:n+1];

            annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics={Line(
                points={{-80,-50},{-80,50},{80,-50},{80,50}},
                color={0,0,255},
                smooth=Smooth.None,
                thickness=1)}));
          end PartialPressureLoss;

          model NominalPressureLoss
        "NominalPressureLoss: Simple pressure loss for nominal values"
            extends PartialPressureLoss;
            import Modelica.Constants.pi;

            parameter Medium.AbsolutePressure dp_nominal
          "Nominal pressure loss";
            parameter Medium.MassFlowRate m_flow_nominal
          "Mass flow rate for dp_nominal";

            parameter Boolean from_dp=true
          " = true, use m_flow = f(dp), otherwise dp = f(m_flow)" 
              annotation (Evaluate=true);
            parameter SI.AbsolutePressure dp_small = 1e-3*dp_nominal
          "Within regularization if |dp| < dp_small" 
              annotation(Dialog(enable=from_dp));
            parameter SI.MassFlowRate m_flow_small = 1e-2*m_flow_nominal
          "Within regularization if |m_flow| < m_flow_small" 
              annotation(Dialog(enable=not from_dp));

            parameter Boolean mixingStreamProperties = true
          "= false to use flow-dependent upstream density and viscosity" 
               annotation(Dialog(group="Advanced"), Evaluate=true);

            parameter Boolean use_d_nominal = false
          "= true, if d_nominal is used, otherwise computed from medium" 
              annotation(Dialog(group="Advanced"), Evaluate=true);
            parameter SI.Density d_nominal = Medium.density_pTX(Medium.p_default, Medium.T_default, Medium.X_default)
          "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
              annotation(Dialog(group="Advanced", enable=use_d_nominal));

            parameter Boolean use_eta_nominal = false
          "= true, if eta_nominal is used, otherwise computed from medium" 
               annotation(Dialog(group="Advanced"), Evaluate=true);
            parameter SI.DynamicViscosity eta_nominal = Medium.dynamicViscosity(
                                                           Medium.setState_pTX(
                                                               Medium.p_default, Medium.T_default, Medium.X_default))
          "Nominal dynamic viscosity (only for m_flow_turbulent and dp_turbulent)"
                  annotation(Dialog(group="Advanced"), Dialog(enable=use_eta_nominal));

            final parameter Boolean constantPressureLossCoefficient=
               use_d_nominal
          "= true if the pressure loss does not depend on fluid states" 
               annotation(Evaluate=true);

            final parameter Boolean continuousFlowReversal=
               mixingStreamProperties
               or constantPressureLossCoefficient
               or not allowFlowReversal
          "= true if the pressure loss is continuous around zero flow" 
               annotation(Evaluate=true);

            final parameter SI.Area[n+1] crossArea_h = {if crossArea[i] > 0 then crossArea[i] else pi/4*2.54e-2^2 for i in 1:n+1}
          "used perimeter, default: diameter of one inch";
            final parameter SI.Length[n] perimeter_h = {if perimeter[i] > 0 then perimeter[i] else pi*2.54e-2 for i in 1:n}
          "used perimeter, default: diameter of one inch";
            final parameter SI.Length[n] roughness_h = {if roughness[i] > 0 then roughness[i] else 2.5e-5 for i in 1:n}
          "used roughness, default: smooth steel pipe";

            SI.Density[n+1] d = if use_d_nominal then fill(d_nominal, n+1) else Medium.density(state);
            SI.Density[n] d_act "actual density per flow segment";

            SI.DynamicViscosity[n+1] eta = if use_eta_nominal then fill(eta_nominal, n+1) else Medium.dynamicViscosity(state);
            SI.DynamicViscosity[n] eta_act "actual viscosity per flow segment";

            SI.Length[n] diameter_h "hydraulic diameter for nominal values";
            SI.Length[n] length_nominal
          "length resulting from nominal pressure loss and geometry";
            Real[n] k_inv "coefficient for quadratic flow";
            Real[n] zeta "coefficient for quadratic flow";

            // Reynolds
            parameter Boolean show_Re = false
          "= true, if Reynolds number, m_flow_turbulent and dp_turbulent are included for plotting"
               annotation (Evaluate=true, Dialog(group="Advanced"));
            SI.ReynoldsNumber[n] Re=Modelica_Fluid.Utilities.ReynoldsNumber_m_flow(
                m_flow/nParallel,
                eta_act,
                diameter_h) if show_Re "Reynolds numbers of pipe flow";
            constant Real Re_turbulent = 4000 "Start of turbulent regime";
            Medium.MassFlowRate[n] m_flow_turbulent=
                {nParallel*(pi/4)*diameter_h[i]*eta_act[i]*Re_turbulent for i in 1:n} if 
                   show_Re "Start of turbulent flow";
            Medium.AbsolutePressure[n] dp_turbulent=
                {(eta_act[i]*diameter_h[i]*pi/4)^2*Re_turbulent^2/(k_inv[i]*d_act[i]) for i in 1:n} if 
                   show_Re "Start of turbulent flow";

          equation
            if not allowFlowReversal then
              d_act = d[1:n];
              eta_act = eta[1:n];
            elseif mixingStreamProperties then
              d_act = 0.5*(d[1:n] + d[2:n+1]);
              eta_act = 0.5*(eta[1:n] + eta[2:n+1]);
            else
              // these values are only used to obtain length_nominal;
              // this is why no events are raised ...
              // the regularization might need to be extended
              if from_dp then
                for i in 1:n loop
                  d_act[i] = noEvent(if dp[i] > 0 then d[i] else d[i+1]);
                  eta_act[i] = noEvent(if dp[i] > 0 then eta[i] else eta[i+1]);
                end for;
              else
                for i in 1:n loop
                  d_act[i] = noEvent(if m_flow[i] > 0 then d[i] else d[i+1]);
                  eta_act[i] = noEvent(if m_flow[i] > 0 then eta[i] else eta[i+1]);
                end for;
              end if;
            end if;
            if continuousFlowReversal and 
               ((from_dp and dp_small >= dp_nominal and dp_nominal > 0)
                 or (not from_dp and m_flow_small >= m_flow_nominal and m_flow_nominal > 0)) then
              // simple linear (laminar) flow
              dp = g*height_ab/n*d_act + dp_nominal/m_flow_nominal*m_flow*nParallel;
            elseif continuousFlowReversal then
              // simple regularization for laminar flow
              if from_dp then
                m_flow = Modelica_Fluid.Pipes.BaseClasses.WallFriction.QuadraticTurbulent.massFlowRate_dp(
                  dp - g*height_ab/n*d_act,
                  d_act,
                  d_act,
                  eta_act,
                  eta_act,
                  length_nominal,
                  diameter_h,
                  roughness_h,
                  dp_small)*nParallel;
              else
                dp = Modelica_Fluid.Pipes.BaseClasses.WallFriction.QuadraticTurbulent.pressureLoss_m_flow(
                  m_flow/nParallel,
                  d_act,
                  d_act,
                  eta_act,
                  eta_act,
                  length_nominal,
                  diameter_h,
                  roughness_h,
                  m_flow_small/nParallel) + g*height_ab/n*d_act;
              end if;
            else
              // regularization for discontinuous flow reversal and static head
              if from_dp then
                m_flow = Modelica_Fluid.Pipes.BaseClasses.WallFriction.QuadraticTurbulent.massFlowRate_dp_staticHead(
                  dp,
                  d[1:n],
                  d[2:n+1],
                  eta[1:n],
                  eta[2:n+1],
                  length_nominal,
                  diameter_h,
                  g*height_ab/n,
                  roughness_h,
                  dp_small/n)*nParallel;
               else
                dp = Modelica_Fluid.Pipes.BaseClasses.WallFriction.QuadraticTurbulent.pressureLoss_m_flow_staticHead(
                  m_flow/nParallel,
                  d[1:n],
                  d[2:n+1],
                  eta[1:n],
                  eta[2:n+1],
                  length_nominal,
                  diameter_h,
                  g*height_ab/n,
                  roughness_h,
                  m_flow_small/nParallel);
               end if;
            end if;
            // Inverse parameterization for WallFriction.QuadraticTurbulent
            // Note: the code should be shared with the WallFriction.QuadraticTurbulent model,
            //       but this required a re-design of the WallFriction interfaces ...
            //   zeta = (length_nominal/diameter)/(2*Math.log10(3.7 /(roughness/diameter)))^2;
            //   k_inv = (pi*diameter*diameter)^2/(8*zeta);
            //   k = d*k_inv "Factor in m_flow = sqrt(k*dp)";
            // and for WallFriction.Laminar.massFlowRate_dp
            //   m_flow = dp*pi*diameter_h^4*d_act/(128*length_laminar*eta_act);
            for i in 1:n loop
              diameter_h[i] = 4*(crossArea_h[i]+crossArea_h[i+1])/2/perimeter_h[i];
              k_inv[i] = (m_flow_nominal/nParallel)^2/((dp_nominal-g*height_ab*d_act[i])/n)/d_act[i];
              zeta[i] = (pi*diameter_h[i]*diameter_h[i])^2/(8*k_inv[i]);
              length_nominal[i] =
                zeta[i]*diameter_h[i]*(2*Modelica.Math.log10(3.7 /(roughness_h[i]/diameter_h[i])))^2;
            end for;

            annotation (Documentation(info="<html>
<p>
This model defines the pressure loss assuming turbulent or laminar flow for 
specified <tt>dp_nominal</tt> and <tt>m_flow_nominal</tt>. 
It takes into account the fluid density of each flow segment and 
obtaines appropriate <tt>length_nominal</tt> values   
for an inverse parameterization of the 
<a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PressureLoss.QuadraticTurbulentFlow\">
          QuadraticTurbulentFlow</a>
model. Per default the upstream and downstream densities are averaged with the setting <tt>mixingStreamProperties = true</tt>,
in order to avoid discontinuous <tt>length_nominal</tt> values in the case of flow reversal.
</p>
<p>
The parameters <tt>dp_small</tt> and <tt>m_flow_small</tt> can be adjusted to account for the laminar zone 
and to numerically regularize the model around zero flow. 
The simplest linear (laminar) pressure loss correlation is obtained with <tt>dp_small >= dp_nominal</tt> for
<tt>from_dp = true</tt> and a continuous flow reversal 
(<tt>mixingStreamProperties = true</tt> or <tt>use_d_nominal = true</tt> or <tt>allowFlowReversal = false</tt>).
</p>
<p>
The geometry parameters <tt>crossArea</tt>, <tt>perimeter</tt> and <tt>roughness</tt> are taken into account if specified. 
Otherwise a smooth steel pipe with a diameter of one inch is used as default for the definition of <tt>length_nominal</tt>. 
The geometry does not effect simulation results of this nominal pressure loss model.
As the geometry is specified however, the internally calculated <tt>m_flow_turbulent</tt> and <tt>dp_turbulent</tt> 
become meaningful and can be related to <tt>m_flow_small</tt> and <tt>dp_small</tt>. 
</p>
<p>
<b>Optional Variables if show_Re</b>
</p>
<table border=1 cellspacing=0 cellpadding=2>
<tr><th><b>Type</b></th><th><b>Name</b></th><th><b>Description</b></th></tr>
<tr><td>ReynoldsNumber</td><td>Re[n]</td>
    <td>Reynolds numbers of pipe flow per flow segment</td></tr> 
<tr><td>MassFlowRate</td><td>m_flow_turbulent[n]</td>
    <td>mass flow rates at start of turbulent region (Re_turbulent=4000)</td></tr>
<tr><td>AbsolutePressure</td><td>dp_turbulent[n]</td>
    <td>pressure losses corresponding to m_flow_turbulent</td></tr>
</table>
</html>", revisions="<html>
<ul>
<li><i>6 Dec 2008</i>
    by Ruediger Franke</a>:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
          end NominalPressureLoss;

          model QuadraticTurbulentFlow
        "Pipe wall friction in the quadratic turbulent regime (simple characteristic, eta not used)"
           extends WallFrictionPressureLoss(
              redeclare package WallFriction = 
              Modelica_Fluid.Pipes.BaseClasses.WallFriction.QuadraticTurbulent);
            annotation (Documentation(info="<html>
<p>
This component defines only the quadratic turbulent regime of wall friction:
dp = k*m_flow*|m_flow|, where \"k\" depends on density and the roughness
of the pipe and is not a function of the Reynolds number.
This relationship is only valid for large Reynolds numbers.
</p>
 
<p>
In the following figure the complete friction regime is shown.
This component describes only the asymptotic behaviour for large
Reynolds numbers, i.e., the values at the right ordinate where
&lambda; is constant.
</p>
 
<img src=\"../Images/Components/PipeFriction1.png\">
 
</html>"));
          end QuadraticTurbulentFlow;

          model WallFrictionPressureLoss
        "WallFrictionPressureLoss: Pipe wall friction in the whole regime (default: detailed characteristic)"
            extends PartialPressureLoss;

            replaceable package WallFriction = 
              Modelica_Fluid.Pipes.BaseClasses.WallFriction.Detailed 
                constrainedby
          Modelica_Fluid.Pipes.BaseClasses.WallFriction.PartialWallFriction
          "Wall friction model" 
                annotation(choicesAllMatching=true,editButton=false);

            parameter Medium.AbsolutePressure dp_nominal = 1e3
          "Nominal pressure loss, only to determine dp_small";
            parameter Medium.MassFlowRate m_flow_nominal = 1
          "Mass flow rate for dp_nominal, only to determine m_flow_small";
            parameter Boolean from_dp=true
          " = true, use m_flow = f(dp), otherwise dp = f(m_flow)" 
              annotation (Evaluate=true);
            parameter SI.AbsolutePressure dp_small = 1e-3*dp_nominal
          "Within regularization if |dp| < dp_small (may be wider for large discontinuities in static head)"
              annotation(Dialog(enable=from_dp and WallFriction.use_dp_small));
            parameter SI.MassFlowRate m_flow_small = 1e-2*m_flow_nominal
          "Within regularization if |m_flow| < m_flow_small (may be wider for large discontinuities in static head)"
              annotation(Dialog(enable=not from_dp and WallFriction.use_m_flow_small));

            parameter Boolean mixingStreamProperties = false
          "= true to use average density and viscosity across flow segments" 
               annotation(Dialog(group="Advanced"), Evaluate=true);

            parameter Boolean use_d_nominal = false
          "= true, if d_nominal is used, otherwise computed from medium" 
               annotation(Dialog(group="Advanced"), Evaluate=true);
            parameter SI.Density d_nominal = Medium.density_pTX(Medium.p_default, Medium.T_default, Medium.X_default)
          "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
              annotation(Dialog(group="Advanced", enable=use_d_nominal));

            parameter Boolean use_eta_nominal = false
          "= true, if eta_nominal is used, otherwise computed from medium" 
               annotation(Dialog(group="Advanced"), Evaluate=true);
            parameter SI.DynamicViscosity eta_nominal = Medium.dynamicViscosity(
                                                           Medium.setState_pTX(
                                                               Medium.p_default, Medium.T_default, Medium.X_default))
          "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
              annotation(Dialog(group="Advanced", enable=use_eta_nominal));

            final parameter Boolean constantPressureLossCoefficient=
               use_d_nominal and (use_eta_nominal or not WallFriction.use_eta)
          "= true if the pressure loss does not depend on fluid states" 
               annotation(Evaluate=true);
            final parameter Boolean continuousFlowReversal=
               mixingStreamProperties
               or constantPressureLossCoefficient
               or not allowFlowReversal
          "= true if the pressure loss is continuous around zero flow" 
               annotation(Evaluate=true);

            parameter Boolean show_Re = false
          "= true, if Reynolds number is included for plotting" 
               annotation (Evaluate=true, Dialog(group="Advanced"));

            SI.ReynoldsNumber[n] Re=Modelica_Fluid.Utilities.ReynoldsNumber_m_flow(
                m_flow/nParallel,
                eta_act,
                diameter) if show_Re "Reynolds numbers of pipe flow";

            // internal variables
            SI.Diameter[n] diameter = {4*(crossArea[i]+crossArea[i+1])/2/perimeter[i] for i in 1:n}
          "Hydraulic diameter";
            SI.Density[n+1] d = if use_d_nominal then fill(d_nominal, n+1) else Medium.density(state);
            SI.Density[n] d_act "Actual density per segment";
            SI.DynamicViscosity[n+1] eta = if not WallFriction.use_eta then fill(1e-10, n+1) else 
                                        (if use_eta_nominal then fill(eta_nominal, n+1) else Medium.dynamicViscosity(state));
            SI.DynamicViscosity[n] eta_act "Actual viscosity per segment";

          equation
            if not allowFlowReversal or constantPressureLossCoefficient then
              d_act = d[1:n];
              eta_act = eta[1:n];
            elseif mixingStreamProperties then
              d_act = 0.5*(d[1:n] + d[2:n+1]);
              eta_act = 0.5*(eta[1:n] + eta[2:n+1]);
            elseif show_Re then
              if from_dp then
                for i in 1:n loop
                  d_act[i] = noEvent(if dp[i] > 0 then d[i] else d[i+1]);
                  eta_act[i] = noEvent(if dp[i] > 0 then eta[i] else eta[i+1]);
                end for;
              else
                for i in 1:n loop
                  d_act[i] = noEvent(if m_flow[i] > 0 then d[i] else d[i+1]);
                  eta_act[i] = noEvent(if m_flow[i] > 0 then eta[i] else eta[i+1]);
                end for;
              end if;
            else // not used for detailed regularization
              d_act = zeros(n);
              eta_act = zeros(n);
            end if;
            if continuousFlowReversal then
              // simple regularization for laminar flow
              if from_dp and not WallFriction.dp_is_zero then
                m_flow = WallFriction.massFlowRate_dp(
                  dp - g*height_ab/n*d_act,
                  d_act,
                  d_act,
                  eta_act,
                  eta_act,
                  length,
                  diameter,
                  roughness,
                  dp_small)*nParallel;
              else
                dp = WallFriction.pressureLoss_m_flow(
                  m_flow/nParallel,
                  d_act,
                  d_act,
                  eta_act,
                  eta_act,
                  length,
                  diameter,
                  roughness,
                  m_flow_small/nParallel) + g*height_ab/n*d_act;
              end if;
            else
              // regularization for discontinuous flow reversal and static head
              if from_dp and not WallFriction.dp_is_zero then
                m_flow = WallFriction.massFlowRate_dp_staticHead(
                  dp,
                  d[1:n],
                  d[2:n+1],
                  eta[1:n],
                  eta[2:n+1],
                  length,
                  diameter,
                  g*height_ab/n,
                  roughness,
                  dp_small/n)*nParallel;
              else
                dp = WallFriction.pressureLoss_m_flow_staticHead(
                  m_flow/nParallel,
                  d[1:n],
                  d[2:n+1],
                  eta[1:n],
                  eta[2:n+1],
                  length,
                  diameter,
                  g*height_ab/n,
                  roughness,
                  m_flow_small/nParallel);
              end if;
            end if;

              annotation (Documentation(info="<html>
<p>
This model describes pressure losses due to <b>wall friction</b> in a pipe
and due to gravity.
It is assumed that no mass or energy is stored in the pipe. 
Correlations of different complexity and validity can be
seleted via the replaceable package <b>WallFriction</b> (see parameter menu below).
The details of the pipe wall friction model are described in the
<a href=\"Modelica://Modelica_Fluid.UsersGuide.ComponentDefinition.WallFriction\">UsersGuide</a>.
Basically, different variants of the equation
</p>
 
<pre>
   dp = &lambda;(Re,<font face=\"Symbol\">D</font>)*(L/D)*&rho;*v*|v|/2
</pre>
 
<p>
are used, where the friction loss factor &lambda; is shown
in the next figure:
</p>
 
<img src=\"../Images/Components/PipeFriction1.png\">
 
<p>
By default, the correlations are computed with media data
at the actual time instant.
In order to reduce non-linear equation systems, the parameters
<b>use_eta_nominal</b> and <b>use_d_nominal</b> provide the option
to compute the correlations with constant media values
at the desired operating point. This might speed-up the
simulation and/or might give a more robust simulation.
</p>
<p>
<b>Optional Variables if show_Re</b>
</p>
<table border=1 cellspacing=0 cellpadding=2>
<tr><th><b>Type</b></th><th><b>Name</b></th><th><b>Description</b></th></tr>
<tr><td>ReynoldsNumber</td><td>Re[n]</td>
    <td>Reynolds numbers of pipe flow per flow segment</td></tr> 
<tr><td>MassFlowRate</td><td>m_flow_turbulent[n]</td>
    <td>mass flow rates at start of turbulent region (Re_turbulent=4000)</td></tr>
</table>
</html>"),    Diagram(coordinateSystem(
                  preserveAspectRatio=false,
                  extent={{-100,-100},{100,100}},
                  grid={1,1}), graphics={
              Rectangle(
                extent={{-100,64},{100,-64}},
                lineColor={0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Backward),
              Rectangle(
                extent={{-100,50},{100,-49}},
                lineColor={0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Line(
                points={{-60,-49},{-60,50}},
                color={0,0,255},
                arrow={Arrow.Filled,Arrow.Filled}),
              Text(
                extent={{-50,16},{6,-10}},
                lineColor={0,0,255},
                fillColor={0,0,255},
                fillPattern=FillPattern.Solid,
                textString="diameter"),
              Line(
                points={{-100,74},{100,74}},
                color={0,0,255},
                arrow={Arrow.Filled,Arrow.Filled}),
              Text(
                extent={{-34,92},{34,74}},
                lineColor={0,0,255},
                fillColor={0,0,255},
                fillPattern=FillPattern.Solid,
                textString="length")}));
          end WallFrictionPressureLoss;

    end PressureLoss;

  package HeatTransfer
    partial model PartialHeatTransfer
        "base class for any pipe heat transfer correlation"

      // Parameters
      replaceable package Medium=Modelica.Media.Interfaces.PartialMedium 
        annotation(Dialog(tab="Internal Interface", enable=false));
      parameter Integer n=1 "Number of heat transfer segments" 
        annotation(Dialog(tab="Internal Interface", enable=false), Evaluate=true);
      parameter Real nParallel "Number of parallel pipes" 
        annotation(Dialog(tab="Internal Interface", enable=false));
      parameter SI.Length[n] length "Length of segments along flow path" 
        annotation(Dialog(tab="Internal Interface", enable=false));
      parameter SI.Area[n+1] crossArea "Cross flow area at segment boundaries" 
        annotation(Dialog(tab="Internal Interface", enable=false));
      parameter SI.Length[n] perimeter
          "Mean perimeter for heat transfer area and hydraulic diameter" 
        annotation(Dialog(tab="Internal Interface", enable=false));
      parameter SI.Height[n] roughness(each min=0)
          "Average height of surface asperities" 
          annotation(Dialog(tab="Internal Interface", enable=false));
      final parameter SI.Length[n] diameter = {4*(crossArea[i]+crossArea[i+1])/2/perimeter[i] for i in 1:n}
          "Hydraulic diameter";
      final parameter SI.Area[n] area = {perimeter[i]*length[i] for i in 1:n}
          "Heat transfer area";
      parameter Boolean use_fluidHeatPort=false
          "= true to use fluidHeatPort instead of output Q_flow" 
        annotation(Dialog(tab="Internal Interface", enable=false));

      // Inputs provided to heat transfer model
      input Medium.ThermodynamicState[n] state;
      input SI.MassFlowRate[n] m_flow;

      // Output defined by heat transfer model
      output SI.HeatFlowRate[n] Q_flow "Heat flow rates per tube";

      // Heat ports
      Modelica_Fluid.Interfaces.HeatPorts_a[n] wallHeatPort "Heat port to wall"
        annotation (Placement(transformation(extent={{-10,60},{10,80}},
                rotation=0), iconTransformation(extent={{-20,60},{20,80}})));

      Modelica_Fluid.Interfaces.HeatPorts_b[n] fluidHeatPort if use_fluidHeatPort
          "Thermal port to fluid" 
        annotation (Placement(transformation(extent={{-10,-70},{10,-50}},
                rotation=0), iconTransformation(extent={{-20,-70},{20,-50}})));

      // Internal variables
      SI.Temperature[n] T;
      Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow[n]
          prescribedHeatFlow "Needed to connect to conditional connector" 
        annotation (Placement(transformation(extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={0,-10})));

    equation
      T = Medium.temperature(state);
      wallHeatPort.Q_flow = Q_flow;
      if use_fluidHeatPort then
        prescribedHeatFlow.Q_flow = Q_flow;
      else
        prescribedHeatFlow.port.T = T;
      end if;

      connect(prescribedHeatFlow.port, fluidHeatPort) annotation (Line(
          points={{-1.83697e-015,-20},{0,-20},{0,-60}},
          color={191,0,0},
          smooth=Smooth.None));

      annotation (Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
                  -100},{100,100}}), graphics={Ellipse(
                extent={{-60,64},{60,-56}},
                lineColor={0,0,0},
                fillPattern=FillPattern.Sphere,
                fillColor={232,0,0}), Text(
                extent={{-38,26},{40,-14}},
                lineColor={0,0,0},
                fillPattern=FillPattern.Sphere,
                fillColor={232,0,0},
                textString="%name")}),
                              Documentation(info="<html>
Base class for heat transfer models that can be used in distributed pipe models.
</html>"),
        Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                  {100,100}}),
                        graphics));
    end PartialHeatTransfer;

    partial model PartialPipeFlowHeatTransfer
        "Base class for pipe heat transfer correlation in terms of Nusselt numberheat transfer in a circular pipe for laminar and turbulent one-phase flow"
      extends PartialHeatTransfer;
      parameter SI.CoefficientOfHeatTransfer alpha0=100
          "guess value for heat transfer coefficients";
      SI.CoefficientOfHeatTransfer[n] alpha(each start=alpha0)
          "CoefficientOfHeatTransfer";
      Real[n] Re "Reynolds number";
      Real[n] Pr "Prandtl number";
      Real[n] Nu "Nusselt number";
      SI.DynamicViscosity[n] eta "Dynamic viscosity";
      SI.ThermalConductivity[n] lambda "Thermal conductivity";
    equation
      eta=Medium.dynamicViscosity(state);
      lambda=Medium.thermalConductivity(state);
      Pr = Medium.prandtlNumber(state);
      Re = CharacteristicNumbers.ReynoldsNumber(m_flow/nParallel, diameter, (crossArea[1:n]+crossArea[2:n+1])/2, eta);
      Nu = CharacteristicNumbers.NusseltNumber(alpha, diameter, lambda);
      Q_flow={alpha[i]*area[i]*(wallHeatPort[i].T - T[i])*nParallel for i in 1:n};
        annotation (Documentation(info="<html>
Base class for heat transfer models that are expressed in terms of the Nusselt number and which can be used in distributed pipe models.
</html>"));
    end PartialPipeFlowHeatTransfer;

    model NoHeatTransfer
        "NoHeatTransfer: No heat transfer assuming perfect isolation"
      extends PartialHeatTransfer;
    equation
      Q_flow = zeros(n);
      annotation(Documentation(info="<html>
Ideal heat transfer without thermal resistance.
</html>"));
    end NoHeatTransfer;

    model IdealHeatTransfer
        "IdealHeatTransfer: Ideal heat transfer without thermal resistance"
      extends PartialHeatTransfer;
    equation
      T = wallHeatPort.T;
      annotation(Documentation(info="<html>
Ideal heat transfer without thermal resistance.
</html>"));
    end IdealHeatTransfer;

    model ConstantHeatTransfer
        "ConstantHeatTransfer: Constant heat transfer coefficient"
      extends PartialHeatTransfer;
      parameter SI.CoefficientOfHeatTransfer alpha0=200
          "heat transfer coefficient";
      annotation(Documentation(info="<html>
Simple heat transfer correlation with constant heat transfer coefficient, used as default component in <a distributed pipe models.
</html>"));
    equation
      Q_flow = {alpha0*area[i]*(wallHeatPort[i].T - T[i])*nParallel for i in 1:n};
    end ConstantHeatTransfer;
    annotation (Documentation(info="<html>
Heat transfer correlations for pipe models
</html>"));

    model LocalPipeFlowHeatTransfer
        "LocalPipeFlowHeatTransfer: Laminar and turbulent forced convection in pipes, local coefficients"
      extends PartialPipeFlowHeatTransfer;
      protected
      Real[n] Nu_turb "Nusselt number for turbulent flow";
      Real[n] Nu_lam "Nusselt number for laminar flow";
      Real Nu_1;
      Real[n] Nu_2;
      Real[n] Xi;
    equation
      Nu_1=3.66;
      for i in 1:n loop
       Nu_turb[i]=smooth(0,(Xi[i]/8)*abs(Re[i])*Pr[i]/(1+12.7*(Xi[i]/8)^0.5*(Pr[i]^(2/3)-1))*(1+1/3*(diameter[i]/length[i]/(if m_flow[i]>=0 then (i-0.5) else (n-i+0.5)))^(2/3)));
       Xi[i]=(1.8*Modelica.Math.log10(max(1e-10,Re[i]))-1.5)^(-2);
       Nu_lam[i]=(Nu_1^3+0.7^3+(Nu_2[i]-0.7)^3)^(1/3);
       Nu_2[i]=smooth(0,1.077*(abs(Re[i])*Pr[i]*diameter[i]/length[i]/(if m_flow[i]>=0 then (i-0.5) else (n-i+0.5)))^(1/3));
       Nu[i]=spliceFunction(Nu_turb[i], Nu_lam[i], Re[i]-6150, 3850);
      end for;
      annotation (Documentation(info="<html>
Heat transfer model for laminar and turbulent flow in pipes. Range of validity:
<ul>
<li>fully developed pipe flow</li>
<li>forced convection</li>
<li>one phase Newtonian fluid</li>
<li>(spatial) constant wall temperature in the laminar region</li>
<li>0 &le; Re &le; 1e6, 0.6 &le; Pr &le; 100, d/L &le; 1</li>
<li>The correlation holds for non-circular pipes only in the turbulent region. Use diameter=4*crossArea/perimeter as characteristic length.</li>
</ul>
The correlation takes into account the spatial position along the pipe flow, which changes discontinuously at flow reversal. However, the heat transfer coefficient itself is continuous around zero flow rate, but not its derivative.
<h4><font color=\"#008000\">References</font></h4>
 
<dl><dt>Verein Deutscher Ingenieure (1997):</dt>
    <dd><b>VDI W&auml;rmeatlas</b>.
         Springer Verlag, Ed. 8, 1997.</dd>
</dl>
</html>"));
    end LocalPipeFlowHeatTransfer;
  end HeatTransfer;

    package CharacteristicNumbers
      function ReynoldsNumber
        input SI.MassFlowRate m_flow "Mass flow rate";
        input SI.Length d_ch "Characteristic length (hyd. diam. in pipes)";
        input SI.Area A "Cross sectional area";
        input SI.DynamicViscosity eta "Dynamic viscosity";
        output SI.ReynoldsNumber Re "Reynolds number";
        annotation (Documentation(info="Calculate Re-Number; Re = mdot*Dhyd/A/eta"),
             Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}),
                  graphics));
      algorithm
        Re := abs(m_flow)*d_ch/A/eta;
      end ReynoldsNumber;

      function NusseltNumber
        input SI.CoefficientOfHeatTransfer alpha "Coefficient of heat transfer";
        input SI.Length d_ch "Characteristic length";
        input SI.ThermalConductivity lambda "Thermal conductivity";
        output SI.NusseltNumber Nu "Nusselt number";
        annotation (Documentation(info="Nusselt number Nu = alpha*d_ch/lambda"));
      algorithm
        Nu := alpha*d_ch/lambda;
      end NusseltNumber;
    end CharacteristicNumbers;

    package WallFriction
      "Different variants for pressure drops due to pipe wall friction"
      partial package PartialWallFriction
        "Partial wall friction characteristic (base package of all wall friction characteristics)"

        annotation (Documentation(info="<html>
 
</html>"));

      // Constants to be set in subpackages
        constant Boolean use_eta = true
          "= true, if eta_a/eta_b are used in function, otherwise value is not used";
        constant Boolean use_roughness = true
          "= true, if roughness is used in function, otherwise value is not used";
        constant Boolean use_dp_small = true
          "= true, if dp_small is used in function, otherwise value is not used";
        constant Boolean use_m_flow_small = true
          "= true, if m_flow_small is used in function, otherwise value is not used";
        constant Boolean dp_is_zero = false
          "= true, if no wall friction is present, i.e., dp = 0 (function massFlowRate_dp() cannot be used)";

      // pressure loss characteristic functions
        replaceable partial function massFlowRate_dp
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction"
          extends Modelica.Icons.Function;

          input SI.Pressure dp "Pressure loss (dp = port_a.p - port_b.p)";
          input SI.Density d_a "Density at port_a";
          input SI.Density d_b "Density at port_b";
          input SI.DynamicViscosity eta_a
            "Dynamic viscosity at port_a (dummy if use_eta = false)";
          input SI.DynamicViscosity eta_b
            "Dynamic viscosity at port_b (dummy if use_eta = false)";
          input SI.Length length "Length of pipe";
          input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
          input SI.Length roughness(min=0) = 2.5e-5
            "Absolute roughness of pipe, with a default for a smooth steel pipe (dummy if use_roughness = false)";
          input SI.AbsolutePressure dp_small = 1
            "Turbulent flow if |dp| >= dp_small (dummy if use_dp_small = false)";

          output SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
        annotation (Documentation(info="<html>
 
</html>"));
        end massFlowRate_dp;

        replaceable partial function massFlowRate_dp_staticHead
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction and static head"
          extends Modelica.Icons.Function;

          input SI.Pressure dp "Pressure loss (dp = port_a.p - port_b.p)";
          input SI.Density d_a "Density at port_a";
          input SI.Density d_b "Density at port_b";
          input SI.DynamicViscosity eta_a
            "Dynamic viscosity at port_a (dummy if use_eta = false)";
          input SI.DynamicViscosity eta_b
            "Dynamic viscosity at port_b (dummy if use_eta = false)";
          input SI.Length length "Length of pipe";
          input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
          input Real g_times_height_ab
            "Gravity times (Height(port_b) - Height(port_a))";
          input SI.Length roughness(min=0) = 2.5e-5
            "Absolute roughness of pipe, with a default for a smooth steel pipe (dummy if use_roughness = false)";
          input SI.AbsolutePressure dp_small=1
            "Turbulent flow if |dp| >= dp_small (dummy if use_dp_small = false)";

          output SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
          annotation (Documentation(info="<html>
 
</html>"));
        end massFlowRate_dp_staticHead;

        replaceable partial function pressureLoss_m_flow
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction"
          extends Modelica.Icons.Function;

          input SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
          input SI.Density d_a "Density at port_a";
          input SI.Density d_b "Density at port_b";
          input SI.DynamicViscosity eta_a
            "Dynamic viscosity at port_a (dummy if use_eta = false)";
          input SI.DynamicViscosity eta_b
            "Dynamic viscosity at port_b (dummy if use_eta = false)";
          input SI.Length length "Length of pipe";
          input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
          input SI.Length roughness(min=0) = 2.5e-5
            "Absolute roughness of pipe, with a default for a smooth steel pipe (dummy if use_roughness = false)";
          input SI.MassFlowRate m_flow_small = 0.01
            "Turbulent flow if |m_flow| >= m_flow_small (dummy if use_m_flow_small = false)";
          output SI.Pressure dp "Pressure loss (dp = port_a.p - port_b.p)";

        annotation (Documentation(info="<html>
 
</html>"));
        end pressureLoss_m_flow;

        replaceable partial function pressureLoss_m_flow_staticHead
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction and static head"
                  extends Modelica.Icons.Function;

          input SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
          input SI.Density d_a "Density at port_a";
          input SI.Density d_b "Density at port_b";
          input SI.DynamicViscosity eta_a
            "Dynamic viscosity at port_a (dummy if use_eta = false)";
          input SI.DynamicViscosity eta_b
            "Dynamic viscosity at port_b (dummy if use_eta = false)";
          input SI.Length length "Length of pipe";
          input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
          input Real g_times_height_ab
            "Gravity times (Height(port_b) - Height(port_a))";
          input SI.Length roughness(min=0) = 2.5e-5
            "Absolute roughness of pipe, with a default for a smooth steel pipe (dummy if use_roughness = false)";
          input SI.MassFlowRate m_flow_small = 0.01
            "Turbulent flow if |m_flow| >= m_flow_small (dummy if use_m_flow_small = false)";
          output SI.Pressure dp "Pressure loss (dp = port_a.p - port_b.p)";

        annotation (Documentation(info="<html>
 
</html>"));
        end pressureLoss_m_flow_staticHead;
      end PartialWallFriction;

      annotation (Documentation(info="<html>
<p>
This package provides functions to compute
pressure losses due to <b>wall friction</b> in a pipe.
Every correlation is defined by a package that is derived
by inheritance from the package WallFriction.PartialWallFriction.
The details of the underlying pipe wall friction model are described in the
<a href=\"Modelica://Modelica_Fluid.UsersGuide.ComponentDefinition.WallFriction\">UsersGuide</a>.
Basically, different variants of the equation
</p>
 
<pre>
   dp = &lambda;(Re,<font face=\"Symbol\">D</font>)*(L/D)*&rho;*v*|v|/2
</pre>
 
<p>
are used, where the friction loss factor &lambda; is shown
in the next figure:
</p>
 
<img src=\"../Images/Components/PipeFriction1.png\">
 
</html>"));
      package NoFriction "No pipe wall friction, no static head"

        annotation (Documentation(info="<html>
<p>
This component sets the pressure loss due to wall friction 
to zero, i.e., it allows to switch off pipe wall friction.
</p>
</html>"));

        extends PartialWallFriction(
                  final use_eta = false,
                  final use_roughness = false,
                  final use_dp_small = false,
                  final use_m_flow_small = false,
                  final dp_is_zero = true);

        redeclare function extends massFlowRate_dp
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction"

          annotation (Documentation(info="<html>
 
</html>"));
        algorithm
          assert(false, "function massFlowRate_dp (option: from_dp=true)
cannot be used for WallFriction.NoFriction. Use instead
function pressureLoss_m_flow (option: from_dp=false)");
        end massFlowRate_dp;

        redeclare function extends pressureLoss_m_flow
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction"

          annotation (Documentation(info="<html>
 
</html>"));
        algorithm
          dp := 0;
        end pressureLoss_m_flow;

        redeclare function extends massFlowRate_dp_staticHead
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction and static head"

          annotation (Documentation(info="<html>
 
</html>"));

        algorithm
          assert(false, "function massFlowRate_dp (option: from_dp=true)
cannot be used for WallFriction.NoFriction. Use instead
function pressureLoss_m_flow (option: from_dp=false)");
        end massFlowRate_dp_staticHead;

        redeclare function extends pressureLoss_m_flow_staticHead
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction and static head"

          annotation (Documentation(info="<html>
 
</html>"));
        /* To include only static head:
protected 
  Real dp_grav_a = g_times_height_ab*d_a 
    "Static head if mass flows in design direction (a to b)";
  Real dp_grav_b = g_times_height_ab*d_b 
    "Static head if mass flows against design direction (b to a)";
*/
        algorithm
        //  dp := Utilities.regStep(m_flow, dp_grav_a, dp_grav_a, m_flow_small);
          dp := 0;
          assert(abs(g_times_height_ab) < Modelica.Constants.small,
           "WallFriction.NoFriction does not consider static head and cannot be used with height_ab<>0!");
        end pressureLoss_m_flow_staticHead;
      end NoFriction;

      package Laminar
        "Pipe wall friction in the laminar regime (linear correlation)"

        annotation (Documentation(info="<html>
<p>
This component defines only the laminar region of wall friction:
dp = k*m_flow, where \"k\" depends on density and dynamic viscosity.
The roughness of the wall does not have an influence on the laminar
flow and therefore argument roughness is ignored.
Since this is a linear relationship, the occuring systems of equations
are usually much simpler (e.g. either linear instead of non-linear).
By using nominal values for density and dynamic viscosity, the 
systems of equations can still further be reduced. 
</p>
 
<p>
In the following figure the complete friction regime is shown.
This component describes only the \"light blue curve\" called
<b>Hagen-Poiseuille</b>.
</p>
 
<img src=\"../Images/Components/PipeFriction1.png\">
 
</html>"));

        extends PartialWallFriction(
                  final use_eta = true,
                  final use_roughness = false,
                  final use_dp_small = false,
                  final use_m_flow_small = false);

        redeclare function extends massFlowRate_dp
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction"

          annotation (Documentation(info="<html>
 
</html>"));
        algorithm
          m_flow :=dp*Modelica.Constants.pi*diameter^4*(d_a + d_b)/(128*length*(eta_a + eta_b));
        end massFlowRate_dp;

        redeclare function extends pressureLoss_m_flow
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction"

          annotation (Documentation(info="<html>
 
</html>"));
        algorithm
          dp := m_flow*128*length*(eta_a + eta_b)/(Modelica.Constants.pi*diameter^4*(d_a + d_b));
        end pressureLoss_m_flow;

        redeclare function extends massFlowRate_dp_staticHead
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction and static head"

          annotation (Documentation(info="<html>
 
</html>"));
        protected
          Real k0inv = Modelica.Constants.pi*diameter^4/(128*length)
            "Constant factor";

          Real dp_grav_a = g_times_height_ab*d_a
            "Static head if mass flows in design direction (a to b)";
          Real dp_grav_b = g_times_height_ab*d_b
            "Static head if mass flows against design direction (b to a)";

          Real dm_flow_ddp_fric_a = k0inv*d_a/eta_a
            "Slope of mass flow rate over dp if flow in design direction (a to b)";
          Real dm_flow_ddp_fric_b = k0inv*d_b/eta_b
            "Slope of mass flow rate over dp if flow against design direction (b to a)";

          Real dp_a=max(dp_grav_a,dp_grav_b)+dp_small
            "Upper end of regularization domain of the m_flow(dp) relation";
          Real dp_b=min(dp_grav_a,dp_grav_b)-dp_small
            "Lower end of regularization domain of the m_flow(dp) relation";

          SI.MassFlowRate m_flow_a
            "Value at upper end of regularization domain";
          SI.MassFlowRate m_flow_b
            "Value at lower end of regularization domain";

          // Properly define zero mass flow conditions
          SI.MassFlowRate m_flow_zero = 0;
          SI.Pressure dp_zero = (dp_grav_a + dp_grav_b)/2;
          Real dm_flow_ddp_fric_zero;
        algorithm
        /*
  dp = 0.5*zeta/(A^2*d) * m_flow * |m_flow|
     = 0.5 * c0/(|m_flow|*(4/pi)/(D_Re*eta)) / ((pi*(D_Re/2)^2)^2*d) * m_flow*|m_flow|
     = 0.5 * c0*(pi/4)*(D_Re*eta) * 16/(pi^2*D_Re^4*d) * m_flow*|m_flow|
     = 2*c0/(pi*D_Re^3) * eta/d * m_flow
     = k0 * eta/d * m_flow
  k0 = 2*c0/(pi*D_Re^3)
*/

          if dp>=dp_a then
            // Positive flow outside regularization
            m_flow := dm_flow_ddp_fric_a*(dp-dp_grav_a);
          elseif dp<=dp_b then
            // Negative flow outside regularization
            m_flow := dm_flow_ddp_fric_b*(dp-dp_grav_b);
          else
            m_flow_a := dm_flow_ddp_fric_a*(dp_a - dp_grav_a);
            m_flow_b := dm_flow_ddp_fric_b*(dp_b - dp_grav_b);

            // Include a properly defined zero mass flow point
            // Obtain a suitable slope from the linear section slope c (value of m_flow is overwritten later)
            (m_flow, dm_flow_ddp_fric_zero) := Utilities.regFun3(dp_zero, dp_b, dp_a, m_flow_b, m_flow_a, dm_flow_ddp_fric_b, dm_flow_ddp_fric_a);
            // Do regularization
            if dp>dp_zero then
              m_flow := Utilities.regFun3(dp, dp_zero, dp_a, m_flow_zero, m_flow_a, dm_flow_ddp_fric_zero, dm_flow_ddp_fric_a);
            else
              m_flow := Utilities.regFun3(dp, dp_b, dp_zero, m_flow_b, m_flow_zero, dm_flow_ddp_fric_b, dm_flow_ddp_fric_zero);
            end if;
          end if;
        /*
  m_flow := if dp<dp_b then dm_flow_ddp_b*(dp-dp_grav_b) else 
              (if dp>dp_a then dm_flow_ddp_a*(dp-dp_grav_a) else 
                Modelica_Fluid.Utilities.regFun3(dp, dp_b, dp_a, dm_flow_ddp_b*(dp_b - dp_grav_b), dm_flow_ddp_a*(dp_a - dp_grav_a), dm_flow_ddp_b, dm_flow_ddp_a));
*/
        end massFlowRate_dp_staticHead;

        redeclare function extends pressureLoss_m_flow_staticHead
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction and static head"

          annotation (Documentation(info="<html>
 
</html>"));
        protected
          Real k0 = 128*length/(Modelica.Constants.pi*diameter^4)
            "Constant factor";

          Real dp_grav_a = g_times_height_ab*d_a
            "Static head if mass flows in design direction (a to b)";
          Real dp_grav_b = g_times_height_ab*d_b
            "Static head if mass flows against design direction (b to a)";

          Real ddp_dm_flow_a = k0*eta_a/d_a
            "Slope of dp over mass flow rate if flow in design direction (a to b)";
          Real ddp_dm_flow_b = k0*eta_b/d_b
            "Slope of dp over mass flow rate if flow against design direction (b to a)";

          Real m_flow_a=if dp_grav_a >= dp_grav_b then m_flow_small else m_flow_small + (dp_grav_b-dp_grav_a)/ddp_dm_flow_a
            "Upper end of regularization domain of the dp(m_flow) relation";
          Real m_flow_b=if dp_grav_a >= dp_grav_b then -m_flow_small else -m_flow_small - (dp_grav_b - dp_grav_a)/ddp_dm_flow_b
            "Lower end of regularization domain of the dp(m_flow) relation";

          SI.Pressure dp_a "Value at upper end of regularization domain";
          SI.Pressure dp_b "Value at lower end of regularization domain";

          // Properly define zero mass flow conditions
          SI.MassFlowRate m_flow_zero = 0;
          SI.Pressure dp_zero = (dp_grav_a + dp_grav_b)/2;
          Real ddp_dm_flow_zero;
        algorithm
        /*
  dp = 0.5*zeta/(A^2*d) * m_flow * |m_flow|
     = 0.5 * c0/(|m_flow|*(4/pi)/(D_Re*eta)) / ((pi*(D_Re/2)^2)^2*d) * m_flow*|m_flow|
     = 0.5 * c0*(pi/4)*(D_Re*eta) * 16/(pi^2*D_Re^4*d) * m_flow*|m_flow|
     = 2*c0/(pi*D_Re^3) * eta/d * m_flow
     = k0 * eta/d * m_flow
  k0 = 2*c0/(pi*D_Re^3)
*/

          if m_flow>=m_flow_a then
            // Positive flow outside regularization
            dp := (ddp_dm_flow_a*m_flow + dp_grav_a);
          elseif m_flow<=m_flow_b then
            // Negative flow outside regularization
            dp := (ddp_dm_flow_b*m_flow + dp_grav_b);
          else
            // Regularization parameters
            dp_a := ddp_dm_flow_a*m_flow_a + dp_grav_a;
            dp_b := ddp_dm_flow_b*m_flow_b + dp_grav_b;
            // Include a properly defined zero mass flow point
            // Obtain a suitable slope from the linear section slope c (value of dp is overwritten later)
            (dp, ddp_dm_flow_zero) := Utilities.regFun3(m_flow_zero, m_flow_b, m_flow_a, dp_b, dp_a, ddp_dm_flow_b, ddp_dm_flow_a);
            // Do regularization
            if m_flow>m_flow_zero then
              dp := Utilities.regFun3(m_flow, m_flow_zero, m_flow_a, dp_zero, dp_a, ddp_dm_flow_zero, ddp_dm_flow_a);
            else
              dp := Utilities.regFun3(m_flow, m_flow_b, m_flow_zero, dp_b, dp_zero, ddp_dm_flow_b, ddp_dm_flow_zero);
            end if;
          end if;
        end pressureLoss_m_flow_staticHead;
      end Laminar;

      package QuadraticTurbulent
        "Pipe wall friction in the quadratic turbulent regime (simple characteristic, eta not used)"

        annotation (Documentation(info="<html>
<p>
This component defines only the quadratic turbulent regime of wall friction:
dp = k*m_flow*|m_flow|, where \"k\" depends on density and the roughness
of the pipe and is no longer a function of the Reynolds number.
This relationship is only valid for large Reynolds numbers.
</p>
 
<p>
In the following figure the complete friction regime is shown.
This component describes only the asymptotic behaviour for large
Reynolds numbers, i.e., the values at the right ordinate where
&lambda; is constant.
</p>
 
<img src=\"../Images/Components/PipeFriction1.png\">
 
</html>"));

        extends PartialWallFriction(
                  final use_eta = false,
                  final use_roughness = true,
                  final use_dp_small = true,
                  final use_m_flow_small = true);

        redeclare function extends massFlowRate_dp
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction"
          import Modelica.Math;
          annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
        protected
          constant Real pi = Modelica.Constants.pi;
          Real zeta;
          Real k_inv;
        algorithm
          /*
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
  */
          assert(roughness > 1.e-10,
                 "roughness > 0 required for quadratic turbulent wall friction characteristic");
          zeta  := (length/diameter)/(2*Math.log10(3.7 /(roughness/diameter)))^2;
          k_inv := (pi*diameter*diameter)^2/(8*zeta);
          m_flow := Modelica_Fluid.Utilities.regRoot2(dp, dp_small, d_a*k_inv, d_b*k_inv);
        end massFlowRate_dp;

        redeclare function extends pressureLoss_m_flow
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction"
          import Modelica.Math;

          annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
        protected
          constant Real pi = Modelica.Constants.pi;
          Real zeta;
          Real k;
        algorithm
          /*
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
  */
          assert(roughness > 1.e-10,
                 "roughness > 0 required for quadratic turbulent wall friction characteristic");
          zeta := (length/diameter)/(2*Math.log10(3.7 /(roughness/diameter)))^2;
          k    := 8*zeta/(pi*diameter*diameter)^2;
          dp   := Modelica_Fluid.Utilities.regSquare2(m_flow, m_flow_small, k/d_a, k/d_b);
        end pressureLoss_m_flow;

        redeclare function extends massFlowRate_dp_staticHead
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction and static head"
          import Modelica.Math;
          annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
        protected
          constant Real pi = Modelica.Constants.pi;
          Real zeta = (length/diameter)/(2*Math.log10(3.7 /(roughness/diameter)))^2;
          Real k_inv = (pi*diameter*diameter)^2/(8*zeta);

          SI.Pressure dp_grav_a = g_times_height_ab*d_a
            "Static head if mass flows in design direction (a to b)";
          SI.Pressure dp_grav_b = g_times_height_ab*d_b
            "Static head if mass flows against design direction (b to a)";

          Real k1 = d_a*k_inv "Factor in m_flow =  sqrt(k1*(dp-dp_grav_a))";
          Real k2 = d_b*k_inv "Factor in m_flow = -sqrt(k2*|dp-dp_grav_b|)";

          Real dp_a=max(dp_grav_a,dp_grav_b)+dp_small
            "Upper end of regularization domain of the m_flow(dp) relation";
          Real dp_b=min(dp_grav_a,dp_grav_b)-dp_small
            "Lower end of regularization domain of the m_flow(dp) relation";

          SI.MassFlowRate m_flow_a
            "Value at upper end of regularization domain";
          SI.MassFlowRate m_flow_b
            "Value at lower end of regularization domain";

          SI.MassFlowRate dm_flow_ddp_fric_a
            "Derivative at upper end of regularization domain";
          SI.MassFlowRate dm_flow_ddp_fric_b
            "Derivative at lower end of regularization domain";

          // Properly define zero mass flow conditions
          SI.MassFlowRate m_flow_zero = 0;
          SI.Pressure dp_zero = (dp_grav_a + dp_grav_b)/2;
          Real dm_flow_ddp_fric_zero;
        algorithm
          /*
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
  */
          assert(roughness > 1.e-10,
                 "roughness > 0 required for quadratic turbulent wall friction characteristic");

          if dp>=dp_a then
            // Positive flow outside regularization
            m_flow := sqrt(k1*(dp-dp_grav_a));
          elseif dp<=dp_b then
            // Negative flow outside regularization
            m_flow := -sqrt(k2*abs(dp-dp_grav_b));
          else
            m_flow_a := sqrt(k1*(dp_a - dp_grav_a));
            m_flow_b := -sqrt(k2*abs(dp_b - dp_grav_b));

            dm_flow_ddp_fric_a := k1/(2*sqrt(k1*(dp_a - dp_grav_a)));
            dm_flow_ddp_fric_b := k2/(2*sqrt(k2*abs(dp_b - dp_grav_b)));
        /*  dm_flow_ddp_fric_a := if abs(dp_a - dp_grav_a)>0 then k1/(2*sqrt(k1*(dp_a - dp_grav_a))) else  Modelica.Constants.inf);
    dm_flow_ddp_fric_b := if abs(dp_b - dp_grav_b)>0 then k2/(2*sqrt(k2*abs(dp_b - dp_grav_b))) else Modelica.Constants.inf; */

            // Include a properly defined zero mass flow point
            // Obtain a suitable slope from the linear section slope c (value of m_flow is overwritten later)
            (m_flow, dm_flow_ddp_fric_zero) := Utilities.regFun3(dp_zero, dp_b, dp_a, m_flow_b, m_flow_a, dm_flow_ddp_fric_b, dm_flow_ddp_fric_a);
            // Do regularization
            if dp>dp_zero then
              m_flow := Utilities.regFun3(dp, dp_zero, dp_a, m_flow_zero, m_flow_a, dm_flow_ddp_fric_zero, dm_flow_ddp_fric_a);
            else
              m_flow := Utilities.regFun3(dp, dp_b, dp_zero, m_flow_b, m_flow_zero, dm_flow_ddp_fric_b, dm_flow_ddp_fric_zero);
            end if;
          end if;
        end massFlowRate_dp_staticHead;

        redeclare function extends pressureLoss_m_flow_staticHead
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction and static head"
          import Modelica.Math;
          annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
        protected
          constant Real pi = Modelica.Constants.pi;
          Real zeta = (length/diameter)/(2*Math.log10(3.7 /(roughness/diameter)))^2;
          Real k = 8*zeta/(pi*diameter*diameter)^2;

          SI.Pressure dp_grav_a = g_times_height_ab*d_a
            "Static head if mass flows in design direction (a to b)";
          SI.Pressure dp_grav_b = g_times_height_ab*d_b
            "Static head if mass flows against design direction (b to a)";

          Real k1 = k/d_a "If m_flow >= 0 then dp = k1*m_flow^2 + dp_grav_a";
          Real k2 = k/d_b "If m_flow < 0 then dp = -k2*m_flow^2 + dp_grav_b";

          Real m_flow_a=if dp_grav_a >= dp_grav_b then m_flow_small else m_flow_small + sqrt((dp_grav_b - dp_grav_a)/k1)
            "Upper end of regularization domain of the dp(m_flow) relation";
          Real m_flow_b=if dp_grav_a >= dp_grav_b then -m_flow_small else -m_flow_small - sqrt((dp_grav_b - dp_grav_a)/k2)
            "Lower end of regularization domain of the dp(m_flow) relation";

          SI.Pressure dp_a "Value at upper end of regularization domain";
          SI.Pressure dp_b "Value at lower end of regularization domain";

          Real ddp_dm_flow_a
            "Derivative of pressure drop with mass flow rate at m_flow_a";
          Real ddp_dm_flow_b
            "Derivative of pressure drop with mass flow rate at m_flow_b";

          // Properly define zero mass flow conditions
          SI.MassFlowRate m_flow_zero = 0;
          SI.Pressure dp_zero = (dp_grav_a + dp_grav_b)/2;
          Real ddp_dm_flow_zero;

        algorithm
          /*
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
  */
          assert(roughness > 1.e-10,
                 "roughness > 0 required for quadratic turbulent wall friction characteristic");

          if m_flow>=m_flow_a then
            // Positive flow outside regularization
            dp := (k1*m_flow^2 + dp_grav_a);
          elseif m_flow<=m_flow_b then
            // Negative flow outside regularization
            dp := (-k2*m_flow^2 + dp_grav_b);
          else
            // Regularization parameters
            dp_a := k1*m_flow_a^2 + dp_grav_a;
            ddp_dm_flow_a := 2*k1*m_flow_a;
            dp_b := -k2*m_flow_b^2 + dp_grav_b;
            ddp_dm_flow_b := -2*k2*m_flow_b;
            // Include a properly defined zero mass flow point
            // Obtain a suitable slope from the linear section slope c (value of dp is overwritten later)
            (dp, ddp_dm_flow_zero) := Utilities.regFun3(m_flow_zero, m_flow_b, m_flow_a, dp_b, dp_a, ddp_dm_flow_b, ddp_dm_flow_a);
            // Do regularization
            if m_flow>m_flow_zero then
              dp := Utilities.regFun3(m_flow, m_flow_zero, m_flow_a, dp_zero, dp_a, ddp_dm_flow_zero, ddp_dm_flow_a);
            else
              dp := Utilities.regFun3(m_flow, m_flow_b, m_flow_zero, dp_b, dp_zero, ddp_dm_flow_b, ddp_dm_flow_zero);
            end if;
          end if;

        end pressureLoss_m_flow_staticHead;
      end QuadraticTurbulent;

      package LaminarAndQuadraticTurbulent
        "Pipe wall friction in the laminar and quadratic turbulent regime (simple characteristic)"

        annotation (Documentation(info="<html>
<p>
This component defines the quadratic turbulent regime of wall friction:
dp = k*m_flow*|m_flow|, where \"k\" depends on density and the roughness
of the pipe and is no longer a function of the Reynolds number.
This relationship is only valid for large Reynolds numbers.
At Re=4000, a polynomial is constructed that approaches
the constant &lambda; (for large Reynolds-numbers) at Re=4000
smoothly and has a derivative at zero mass flow rate that is
identical to laminar wall friction.
</p>
</html>"));

        extends PartialWallFriction(
                  final use_eta = true,
                  final use_roughness = true,
                  final use_dp_small = true,
                  final use_m_flow_small = true);

        import ln = Modelica.Math.log "Logarithm, base e";
        import Modelica.Math.log10 "Logarithm, base 10";
        import Modelica.Math.exp "Exponential function";
        import Modelica.Constants.pi;

        redeclare function extends massFlowRate_dp
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction"
          import Modelica.Math;
          annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
        protected
          constant Real pi=Modelica.Constants.pi;
          constant Real Re_turbulent = 4000 "Start of turbulent regime";
          Real zeta;
          Real k0;
          Real k_inv;
          Real yd0 "Derivative of m_flow=m_flow(dp) at zero";
          SI.AbsolutePressure dp_turbulent;
        algorithm
        /*
Turbulent region:
   Re = m_flow*(4/pi)/(D_Re*eta)  
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
   m_flow_turbulent = (pi/4)*D_Re*eta*Re_turbulent
   dp_turbulent     =  k/d *(D_Re*eta*pi/4)^2 * Re_turbulent^2
 
   The start of the turbulent region is computed with mean values
   of dynamic viscosity eta and density rho. Otherwise, one has
   to introduce different "delta" values for both flow directions.
   In order to simplify the approach, only one delta is used.  
 
Laminar region:
   dp = 0.5*zeta/(A^2*d) * m_flow * |m_flow|
      = 0.5 * c0/(|m_flow|*(4/pi)/(D_Re*eta)) / ((pi*(D_Re/2)^2)^2*d) * m_flow*|m_flow|
      = 0.5 * c0*(pi/4)*(D_Re*eta) * 16/(pi^2*D_Re^4*d) * m_flow*|m_flow|
      = 2*c0/(pi*D_Re^3) * eta/d * m_flow
      = k0 * eta/d * m_flow
   k0 = 2*c0/(pi*D_Re^3)
 
   In order that the derivative of dp=f(m_flow) is continuous 
   at m_flow=0, the mean values of eta and d are used in the
   laminar region: eta/d = (eta_a + eta_b)/(d_a + d_b)
   If data.zetaLaminarKnown = false then eta_a and eta_b are potentially zero
   (because dummy values) and therefore the division is only performed
   if zetaLaminarKnown = true.
*/
          assert(roughness > 1.e-10,
                 "roughness > 0 required for quadratic turbulent wall friction characteristic");
          zeta   := (length/diameter)/(2*Math.log10(3.7 /(roughness/diameter)))^2;
          k0     := 128*length/(pi*diameter^4);
          k_inv  := (pi*diameter*diameter)^2/(8*zeta);
          yd0    := (d_a + d_b)/(k0*(eta_a + eta_b));
          dp_turbulent := ((eta_a + eta_b)*diameter*pi/8)^2*Re_turbulent^2/(k_inv*(d_a+d_b)/2);
          m_flow := Modelica_Fluid.Utilities.regRoot2(dp, dp_turbulent, d_a*k_inv, d_b*k_inv,
                                                      use_yd0=true, yd0=yd0);
        end massFlowRate_dp;

        redeclare function extends pressureLoss_m_flow
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction"
          import Modelica.Math;

          annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
        protected
          constant Real pi=Modelica.Constants.pi;
          constant Real Re_turbulent = 4000 "Start of turbulent regime";
          Real zeta;
          Real k0;
          Real k;
          Real yd0 "Derivative of dp = f(m_flow) at zero";
          SI.MassFlowRate m_flow_turbulent
            "The turbulent region is: |m_flow| >= m_flow_turbulent";

        algorithm
        /*
Turbulent region:
   Re = m_flow*(4/pi)/(D_Re*eta)  
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
   m_flow_turbulent = (pi/4)*D_Re*eta*Re_turbulent
   dp_turbulent     =  k/d *(D_Re*eta*pi/4)^2 * Re_turbulent^2
 
   The start of the turbulent region is computed with mean values
   of dynamic viscosity eta and density rho. Otherwise, one has
   to introduce different "delta" values for both flow directions.
   In order to simplify the approach, only one delta is used.  
 
Laminar region:
   dp = 0.5*zeta/(A^2*d) * m_flow * |m_flow|
      = 0.5 * c0/(|m_flow|*(4/pi)/(D_Re*eta)) / ((pi*(D_Re/2)^2)^2*d) * m_flow*|m_flow|
      = 0.5 * c0*(pi/4)*(D_Re*eta) * 16/(pi^2*D_Re^4*d) * m_flow*|m_flow|
      = 2*c0/(pi*D_Re^3) * eta/d * m_flow
      = k0 * eta/d * m_flow
   k0 = 2*c0/(pi*D_Re^3)
 
   In order that the derivative of dp=f(m_flow) is continuous 
   at m_flow=0, the mean values of eta and d are used in the
   laminar region: eta/d = (eta_a + eta_b)/(d_a + d_b)
*/
          assert(roughness > 1.e-10,
                 "roughness > 0 required for quadratic turbulent wall friction characteristic");
          zeta := (length/diameter)/(2*Math.log10(3.7 /(roughness/diameter)))^2;
          k0   := 128*length/(pi*diameter^4);
          k    := 8*zeta/(pi*diameter*diameter)^2;
          yd0  := k0*(eta_a + eta_b)/(d_a + d_b);
          m_flow_turbulent :=(pi/8)*diameter*(eta_a + eta_b)*Re_turbulent;
          dp :=Modelica_Fluid.Utilities.regSquare2(m_flow, m_flow_turbulent, k/d_a, k/d_b,
                                                   use_yd0=true, yd0=yd0);
        end pressureLoss_m_flow;

        redeclare function extends massFlowRate_dp_staticHead
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction and static head"
          import Modelica.Math;
          annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));

        protected
          Real Delta = roughness/diameter "Relative roughness";
          SI.ReynoldsNumber Re1 = 745*exp(if Delta <= 0.0065 then 1 else 0.0065/Delta)
            "Boundary between laminar regime and transition";
          constant SI.ReynoldsNumber Re2 = 4000
            "Boundary between transition and turbulent regime";

          SI.Pressure dp_a
            "Upper end of regularization domain of the m_flow(dp) relation";
          SI.Pressure dp_b
            "Lower end of regularization domain of the m_flow(dp) relation";

          SI.MassFlowRate m_flow_a
            "Value at upper end of regularization domain";
          SI.MassFlowRate m_flow_b
            "Value at lower end of regularization domain";

          SI.MassFlowRate dm_flow_ddp_fric_a
            "Derivative at upper end of regularization domain";
          SI.MassFlowRate dm_flow_ddp_fric_b
            "Derivative at lower end of regularization domain";

          SI.Pressure dp_grav_a = g_times_height_ab*d_a
            "Static head if mass flows in design direction (a to b)";
          SI.Pressure dp_grav_b = g_times_height_ab*d_b
            "Static head if mass flows against design direction (b to a)";

          // Properly define zero mass flow conditions
          SI.MassFlowRate m_flow_zero = 0;
          SI.Pressure dp_zero = (dp_grav_a + dp_grav_b)/2;
          Real dm_flow_ddp_fric_zero;
        algorithm
          assert(roughness > 1.e-10,
            "roughness > 0 required for quadratic turbulent wall friction characteristic");

          dp_a := max(dp_grav_a, dp_grav_b)+dp_small;
          dp_b := min(dp_grav_a, dp_grav_b)-dp_small;

          if dp>=dp_a then
            // Positive flow outside regularization
            m_flow := Internal.m_flow_of_dp_fric(dp - dp_grav_a, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta);
          elseif dp<=dp_b then
            // Negative flow outside regularization
            m_flow := Internal.m_flow_of_dp_fric(dp-dp_grav_b, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta);
          else
            // Regularization parameters
            (m_flow_a, dm_flow_ddp_fric_a) := Internal.m_flow_of_dp_fric(dp_a-dp_grav_a, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta);
            (m_flow_b, dm_flow_ddp_fric_b) := Internal.m_flow_of_dp_fric(dp_b-dp_grav_b, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta);
            // Include a properly defined zero mass flow point
            // Obtain a suitable slope from the linear section slope c (value of m_flow is overwritten later)
            (m_flow, dm_flow_ddp_fric_zero) := Utilities.regFun3(dp_zero, dp_b, dp_a, m_flow_b, m_flow_a, dm_flow_ddp_fric_b, dm_flow_ddp_fric_a);
            // Do regularization
            if dp>dp_zero then
              m_flow := Utilities.regFun3(dp, dp_zero, dp_a, m_flow_zero, m_flow_a, dm_flow_ddp_fric_zero, dm_flow_ddp_fric_a);
            else
              m_flow := Utilities.regFun3(dp, dp_b, dp_zero, m_flow_b, m_flow_zero, dm_flow_ddp_fric_b, dm_flow_ddp_fric_zero);
            end if;
          end if;
        end massFlowRate_dp_staticHead;

        redeclare function extends pressureLoss_m_flow_staticHead
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction and static head"
          import Modelica.Math;
          annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));

        protected
          Real Delta = roughness/diameter "Relative roughness";
          SI.ReynoldsNumber Re1 = 745*exp(if Delta <= 0.0065 then 1 else 0.0065/Delta)
            "Boundary between laminar regime and transition";
          constant SI.ReynoldsNumber Re2 = 4000
            "Boundary between transition and turbulent regime";

          SI.MassFlowRate m_flow_a
            "Upper end of regularization domain of the dp(m_flow) relation";
          SI.MassFlowRate m_flow_b
            "Lower end of regularization domain of the dp(m_flow) relation";

          SI.Pressure dp_a "Value at upper end of regularization domain";
          SI.Pressure dp_b "Value at lower end of regularization domain";

          SI.Pressure dp_grav_a = g_times_height_ab*d_a
            "Static head if mass flows in design direction (a to b)";
          SI.Pressure dp_grav_b = g_times_height_ab*d_b
            "Static head if mass flows against design direction (b to a)";

          Real ddp_dm_flow_a
            "Derivative of pressure drop with mass flow rate at m_flow_a";
          Real ddp_dm_flow_b
            "Derivative of pressure drop with mass flow rate at m_flow_b";

          // Properly define zero mass flow conditions
          SI.MassFlowRate m_flow_zero = 0;
          SI.Pressure dp_zero = (dp_grav_a + dp_grav_b)/2;
          Real ddp_dm_flow_zero;

        algorithm
          assert(roughness > 1.e-10,
            "roughness > 0 required for quadratic turbulent wall friction characteristic");

          m_flow_a := if dp_grav_a<dp_grav_b then 
            Internal.m_flow_of_dp_fric(dp_grav_b - dp_grav_a, d_a, d_b, eta_a, eta_b, length, diameter,  Re1, Re2, Delta)+m_flow_small else 
            m_flow_small;
          m_flow_b := if dp_grav_a<dp_grav_b then 
            Internal.m_flow_of_dp_fric(dp_grav_a - dp_grav_b, d_a, d_b, eta_a, eta_b, length, diameter,  Re1, Re2, Delta)-m_flow_small else 
            -m_flow_small;

          if m_flow>=m_flow_a then
            // Positive flow outside regularization
            dp := Internal.dp_fric_of_m_flow(m_flow, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta) + dp_grav_a;
          elseif m_flow<=m_flow_b then
            // Negative flow outside regularization
            dp := Internal.dp_fric_of_m_flow(m_flow, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta) + dp_grav_b;
          else
            // Regularization parameters
            (dp_a, ddp_dm_flow_a) := Internal.dp_fric_of_m_flow(m_flow_a, d_a, d_b, eta_a, eta_b, length, diameter,  Re1, Re2, Delta);
            dp_a := dp_a + dp_grav_a "Adding dp_grav to dp_fric to get dp";
            (dp_b, ddp_dm_flow_b) := Internal.dp_fric_of_m_flow(m_flow_b, d_a, d_b, eta_a, eta_b, length, diameter,  Re1, Re2, Delta);
            dp_b := dp_b + dp_grav_b "Adding dp_grav to dp_fric to get dp";
            // Include a properly defined zero mass flow point
            // Obtain a suitable slope from the linear section slope c (value of dp is overwritten later)
            (dp, ddp_dm_flow_zero) := Utilities.regFun3(m_flow_zero, m_flow_b, m_flow_a, dp_b, dp_a, ddp_dm_flow_b, ddp_dm_flow_a);
            // Do regularization
            if m_flow>m_flow_zero then
              dp := Utilities.regFun3(m_flow, m_flow_zero, m_flow_a, dp_zero, dp_a, ddp_dm_flow_zero, ddp_dm_flow_a);
            else
              dp := Utilities.regFun3(m_flow, m_flow_b, m_flow_zero, dp_b, dp_zero, ddp_dm_flow_b, ddp_dm_flow_zero);
            end if;
          end if;
        end pressureLoss_m_flow_staticHead;

        package Internal
          "Functions to calculate mass flow rate from friction pressure drop and vice versa"
          function m_flow_of_dp_fric
            "Calculate mass flow rate as function of pressure drop due to friction"

            input SI.Pressure dp_fric
              "Pressure loss due to friction (dp = port_a.p - port_b.p)";
            input SI.Density d_a "Density at port_a";
            input SI.Density d_b "Density at port_b";
            input SI.DynamicViscosity eta_a
              "Dynamic viscosity at port_a (dummy if use_eta = false)";
            input SI.DynamicViscosity eta_b
              "Dynamic viscosity at port_b (dummy if use_eta = false)";
            input SI.Length length "Length of pipe";
            input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
            input SI.ReynoldsNumber Re1
              "Boundary between laminar regime and transition";
            input SI.ReynoldsNumber Re2
              "Boundary between transition and turbulent regime";
            input Real Delta "Relative roughness";
            output SI.MassFlowRate m_flow
              "Mass flow rate from port_a to port_b";
            output Real dm_flow_ddp_fric
              "Derivative of mass flow rate with dp_fric";
            annotation (smoothOrder=1);
          protected
            SI.DynamicViscosity eta "Upstream viscosity";
            SI.Density d "Upstream density";

            Real zeta;
            Real k0;
            Real k_inv;
            Real dm_flow_ddp_laminar
              "Derivative of m_flow=m_flow(dp) in laminar regime";
            SI.AbsolutePressure dp_fric_turbulent
              "The turbulent region is: |dp_fric| >= dp_fric_turbulent, simple quadratic correlation";
            SI.AbsolutePressure dp_fric_laminar
              "The laminar region is: |dp_fric| <= dp_fric_laminar";
          algorithm
          /*
Turbulent region:
   Re = m_flow*(4/pi)/(D_Re*eta)  
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
   dp_fric_turbulent     =  k/d *(D_Re*eta*pi/4)^2 * Re_turbulent^2  
 
Laminar region:
   dp = 0.5*zeta/(A^2*d) * m_flow * |m_flow|
      = 0.5 * c0/(|m_flow|*(4/pi)/(D_Re*eta)) / ((pi*(D_Re/2)^2)^2*d) * m_flow*|m_flow|
      = 0.5 * c0*(pi/4)*(D_Re*eta) * 16/(pi^2*D_Re^4*d) * m_flow*|m_flow|
      = 2*c0/(pi*D_Re^3) * eta/d * m_flow
      = k0 * eta/d * m_flow
   k0 = 2*c0/(pi*D_Re^3)
*/
            // Determine upstream density and upstream viscosity
            if dp_fric >= 0 then
              d := d_a;
              eta := eta_a;
            else
              d := d_b;
              eta := eta_b;
            end if;
            // Quadratic turbulent
            zeta := (length/diameter)/(2*log10(3.7/(Delta)))^2;
            k_inv := (pi*diameter*diameter)^2/(8*zeta);
            dp_fric_turbulent := sign(dp_fric)*(eta*diameter*pi/4)^2*Re2^2/(k_inv*d);

            // Laminar
            k0 := 128*length/(pi*diameter^4);
            dm_flow_ddp_laminar := d/(k0*eta);
            dp_fric_laminar := sign(dp_fric)*pi*k0*eta^2/d*diameter/4*Re1;

            if abs(dp_fric) > abs(dp_fric_turbulent) then
              m_flow := sign(dp_fric)*sqrt(d*k_inv*abs(dp_fric));
              dm_flow_ddp_fric := 0.5*d*k_inv*(d*k_inv*abs(dp_fric))^(-0.5);
            elseif abs(dp_fric) < abs(dp_fric_laminar) then
              m_flow := dm_flow_ddp_laminar*dp_fric;
              dm_flow_ddp_fric := dm_flow_ddp_laminar;
            else
              // Preliminary testing seems to indicate that the log-log transform is not required here
              (m_flow,dm_flow_ddp_fric) := Utilities.cubicHermite_withDerivative(
                dp_fric, dp_fric_laminar, dp_fric_turbulent, dm_flow_ddp_laminar*dp_fric_laminar,
                sign(dp_fric_turbulent)*sqrt(d*k_inv*abs(dp_fric_turbulent)), dm_flow_ddp_laminar,
                if abs(dp_fric_turbulent)>0 then 0.5*d*k_inv*(d*k_inv*abs(dp_fric_turbulent))^(-0.5) else Modelica.Constants.inf);
            end if;
          end m_flow_of_dp_fric;

          function dp_fric_of_m_flow
            "Calculate pressure drop due to friction as function of mass flow rate"

            input SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
            input SI.Density d_a "Density at port_a";
            input SI.Density d_b "Density at port_b";
            input SI.DynamicViscosity eta_a
              "Dynamic viscosity at port_a (dummy if use_eta = false)";
            input SI.DynamicViscosity eta_b
              "Dynamic viscosity at port_b (dummy if use_eta = false)";
            input SI.Length length "Length of pipe";
            input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
            input SI.ReynoldsNumber Re1
              "Boundary between laminar regime and transition";
            input SI.ReynoldsNumber Re2
              "Boundary between transition and turbulent regime";
            input Real Delta "Relative roughness";
            output SI.Pressure dp_fric
              "Pressure loss due to friction (dp_fric = port_a.p - port_b.p - dp_grav)";
            output Real ddp_fric_dm_flow
              "Derivative of pressure drop with mass flow rate";
            annotation (smoothOrder=1);
          protected
            SI.DynamicViscosity eta "Upstream viscosity";
            SI.Density d "Upstream density";
            Real zeta;
            Real k0;
            Real k;
            Real ddp_fric_dm_flow_laminar
              "Derivative of dp_fric = f(m_flow) at zero";
            SI.MassFlowRate m_flow_turbulent
              "The turbulent region is: |m_flow| >= m_flow_turbulent";
            SI.MassFlowRate m_flow_laminar
              "The laminar region is: |m_flow| <= m_flow_laminar";
          algorithm
          /*
Turbulent region:
   Re = m_flow*(4/pi)/(D_Re*eta)  
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
   m_flow_turbulent = (pi/4)*D_Re*eta*Re_turbulent
 
Laminar region:
   dp = 0.5*zeta/(A^2*d) * m_flow * |m_flow|
      = 0.5 * c0/(|m_flow|*(4/pi)/(D_Re*eta)) / ((pi*(D_Re/2)^2)^2*d) * m_flow*|m_flow|
      = 0.5 * c0*(pi/4)*(D_Re*eta) * 16/(pi^2*D_Re^4*d) * m_flow*|m_flow|
      = 2*c0/(pi*D_Re^3) * eta/d * m_flow
      = k0 * eta/d * m_flow
   k0 = 2*c0/(pi*D_Re^3)
*/
            // Determine upstream density and upstream viscosity
            if m_flow >= 0 then
              d := d_a;
              eta := eta_a;
            else
              d := d_b;
              eta := eta_b;
            end if;

            // Turbulent
            zeta := (length/diameter)/(2*log10(3.7/(Delta)))^2;
            k := 8*zeta/(pi*diameter*diameter)^2;
            m_flow_turbulent := sign(m_flow)*(pi/4)*diameter*eta*Re2;

            // Laminar
            k0 := 128*length/(pi*diameter^4);
            ddp_fric_dm_flow_laminar := k0*eta/d;
            m_flow_laminar := sign(m_flow)*(pi/4)*diameter*eta*Re1;

            if abs(m_flow) > abs(m_flow_turbulent) then
              dp_fric := k/d*m_flow*abs(m_flow);
              ddp_fric_dm_flow := 2*k/d*abs(m_flow);
            elseif abs(m_flow) < abs(m_flow_laminar) then
              dp_fric := ddp_fric_dm_flow_laminar*m_flow;
              ddp_fric_dm_flow := ddp_fric_dm_flow_laminar;
            else
              // Preliminary testing seems to indicate that the log-log transform is not required here
              (dp_fric,ddp_fric_dm_flow) := Utilities.cubicHermite_withDerivative(
                m_flow, m_flow_laminar, m_flow_turbulent, ddp_fric_dm_flow_laminar*m_flow_laminar,
                k/d*m_flow_turbulent*abs(m_flow_turbulent), ddp_fric_dm_flow_laminar, 2*k/d*abs(m_flow_turbulent));
            end if;
          end dp_fric_of_m_flow;
        end Internal;
      end LaminarAndQuadraticTurbulent;

      package Detailed
        "Pipe wall friction in the whole regime (detailed characteristic)"

        annotation (Documentation(info="<html>
<p>
This component defines the complete regime of wall friction.
The details are described in the
<a href=\"Modelica://Modelica_Fluid.UsersGuide.ComponentDefinition.WallFriction\">UsersGuide</a>.
The functional relationship of the friction loss factor &lambda; is
displayed in the next figure. Function massFlowRate_dp() defines the \"red curve\"
(\"Swamee and Jain\"), where as function pressureLoss_m_flow() defines the
\"blue curve\" (\"Colebrook-White\"). The two functions are inverses from 
each other and give slightly different results in the transition region
between Re = 1500 .. 4000, in order to get explicit equations without
solving a non-linear equation.
</p>
 
<img src=\"../Images/Components/PipeFriction1.png\">
 
<p>
Additionally to wall friction, this component properly implements static
head. With respect to the latter, two cases can be distuinguised. In the case
shown next, the change of elevation with the path from a to b has the opposite
sign of the change of density.</p>
 
<img src=\"../Images/Components/PipeFrictionStaticHead_case-a.PNG\">
 
<p>
In the case illustrated second, the change of elevation with the path from a to 
b has the same sign of the change of density.</p>
 
<img src=\"../Images/Components/PipeFrictionStaticHead_case-b.PNG\">
 
 
</html>"));

        extends PartialWallFriction(
                  final use_eta = true,
                  final use_roughness = true,
                  final use_dp_small = true,
                  final use_m_flow_small = true);

        import ln = Modelica.Math.log "Logarithm, base e";
        import Modelica.Math.log10 "Logarithm, base 10";
        import Modelica.Math.exp "Exponential function";
        import Modelica.Constants.pi;

        redeclare function extends massFlowRate_dp
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction"
          import Modelica.Math;
                  annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
        protected
          constant Real pi = Modelica.Constants.pi;
          Real Delta = roughness/diameter "Relative roughness";
          SI.ReynoldsNumber Re1 = (745*Math.exp(if Delta <= 0.0065 then 1 else 0.0065/Delta))^0.97
            "Re leaving laminar curve";
          SI.ReynoldsNumber Re2 = 4000 "Re entering turbulent curve";
          SI.DynamicViscosity eta "Upstream viscosity";
          SI.Density d "Upstream density";
          SI.ReynoldsNumber Re "Reynolds number";
          Real lambda2 "Modified friction coefficient (= lambda*Re^2)";

          function interpolateInRegion2
             input Real Re_turbulent;
             input SI.ReynoldsNumber Re1;
             input SI.ReynoldsNumber Re2;
             input Real Delta;
             input Real lambda2;
             output SI.ReynoldsNumber Re;
             annotation(smoothOrder=1);
            // point lg(lambda2(Re1)) with derivative at lg(Re1)
          protected
            Real x1=Math.log10(64*Re1);
            Real y1=Math.log10(Re1);
            Real yd1=1;

            // Point lg(lambda2(Re2)) with derivative at lg(Re2)
            Real aux1=(0.5/Math.log(10))*5.74*0.9;
            Real aux2=Delta/3.7 + 5.74/Re2^0.9;
            Real aux3=Math.log10(aux2);
            Real L2=0.25*(Re2/aux3)^2;
            Real aux4=2.51/sqrt(L2) + 0.27*Delta;
            Real aux5=-2*sqrt(L2)*Math.log10(aux4);
            Real x2=Math.log10(L2);
            Real y2=Math.log10(aux5);
            Real yd2=0.5 + (2.51/Math.log(10))/(aux5*aux4);

            // Constants: Cubic polynomial between lg(Re1) and lg(Re2)
            Real diff_x=x2 - x1;
            Real m=(y2 - y1)/diff_x;
            Real c2=(3*m - 2*yd1 - yd2)/diff_x;
            Real c3=(yd1 + yd2 - 2*m)/(diff_x*diff_x);
            Real lambda2_1=64*Re1;
            Real dx;
          algorithm
             dx := Math.log10(lambda2/lambda2_1);
             Re := Re1*(lambda2/lambda2_1)^(1 + dx*(c2 + dx*c3));
          end interpolateInRegion2;

        algorithm
          // Determine upstream density, upstream viscosity, and lambda2
          d       := if dp >= 0 then d_a else d_b;
          eta     := if dp >= 0 then eta_a else eta_b;
          lambda2 := abs(dp)*2*diameter^3*d/(length*eta*eta);

          // Determine Re under the assumption of laminar flow
          Re := lambda2/64;

          // Modify Re, if turbulent flow
          if Re > Re1 then
             Re :=-2*sqrt(lambda2)*Math.log10(2.51/sqrt(lambda2) + 0.27*Delta);
             if Re < Re2 then
                Re := interpolateInRegion2(Re, Re1, Re2, Delta, lambda2);
             end if;
          end if;

          // Determine mass flow rate
          m_flow := (pi*diameter/4)*eta*(if dp >= 0 then Re else -Re);
        end massFlowRate_dp;

        redeclare function extends pressureLoss_m_flow
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction"
          import Modelica.Math;
                  annotation (smoothOrder=1, Documentation(info="<html>
 
</html>"));
        protected
          constant Real pi = Modelica.Constants.pi;
          Real Delta = roughness/diameter "Relative roughness";
          SI.ReynoldsNumber Re1 = 745*Math.exp(if Delta <= 0.0065 then 1 else 0.0065/Delta)
            "Re leaving laminar curve";
          SI.ReynoldsNumber Re2 = 4000 "Re entering turbulent curve";
          SI.DynamicViscosity eta "Upstream viscosity";
          SI.Density d "Upstream density";
          SI.ReynoldsNumber Re "Reynolds number";
          Real lambda2 "Modified friction coefficient (= lambda*Re^2)";

          function interpolateInRegion2
             input SI.ReynoldsNumber Re;
             input SI.ReynoldsNumber Re1;
             input SI.ReynoldsNumber Re2;
             input Real Delta;
             output Real lambda2;
             annotation(smoothOrder=1);
            // point lg(lambda2(Re1)) with derivative at lg(Re1)
          protected
            Real x1 = Math.log10(Re1);
            Real y1 = Math.log10(64*Re1);
            Real yd1=1;

            // Point lg(lambda2(Re2)) with derivative at lg(Re2)
            Real aux1=(0.5/Math.log(10))*5.74*0.9;
            Real aux2=Delta/3.7 + 5.74/Re2^0.9;
            Real aux3=Math.log10(aux2);
            Real L2=0.25*(Re2/aux3)^2;
            Real aux4=2.51/sqrt(L2) + 0.27*Delta;
            Real aux5=-2*sqrt(L2)*Math.log10(aux4);
            Real x2 =  Math.log10(Re2);
            Real y2 =  Math.log10(L2);
            Real yd2 = 2 + 4*aux1/(aux2*aux3*(Re2)^0.9);

            // Constants: Cubic polynomial between lg(Re1) and lg(Re2)
            Real diff_x=x2 - x1;
            Real m=(y2 - y1)/diff_x;
            Real c2=(3*m - 2*yd1 - yd2)/diff_x;
            Real c3=(yd1 + yd2 - 2*m)/(diff_x*diff_x);
            Real dx;
          algorithm
             dx := Math.log10(Re/Re1);
             lambda2 := 64*Re1*(Re/Re1)^(1 + dx*(c2 + dx*c3));
          end interpolateInRegion2;
        algorithm
          // Determine upstream density and upstream viscosity
          d       :=if m_flow >= 0 then d_a else d_b;
          eta     :=if m_flow >= 0 then eta_a else eta_b;

          // Determine Re, lambda2 and pressure drop
          Re :=(4/pi)*abs(m_flow)/(diameter*eta);
          lambda2 := if Re <= Re1 then 64*Re else 
                    (if Re >= Re2 then 0.25*(Re/Math.log10(Delta/3.7 + 5.74/Re^0.9))^2 else 
                     interpolateInRegion2(Re, Re1, Re2, Delta));
          dp :=length*eta*eta/(2*d*diameter*diameter*diameter)*
               (if m_flow >= 0 then lambda2 else -lambda2);
        end pressureLoss_m_flow;

        redeclare function extends massFlowRate_dp_staticHead
          "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction and static head"

          annotation (smoothOrder=1);

        protected
          Real Delta = roughness/diameter "Relative roughness";
          SI.ReynoldsNumber Re "Reynolds number";
          SI.ReynoldsNumber Re1 = (745*exp(if Delta <= 0.0065 then 1 else 0.0065/Delta))^0.97
            "Boundary between laminar regime and transition";
          constant SI.ReynoldsNumber Re2 = 4000
            "Boundary between transition and turbulent regime";
          SI.Pressure dp_a
            "Upper end of regularization domain of the m_flow(dp) relation";
          SI.Pressure dp_b
            "Lower end of regularization domain of the m_flow(dp) relation";
          SI.MassFlowRate m_flow_a
            "Value at upper end of regularization domain";
          SI.MassFlowRate m_flow_b
            "Value at lower end of regularization domain";

          SI.MassFlowRate dm_flow_ddp_fric_a
            "Derivative at upper end of regularization domain";
          SI.MassFlowRate dm_flow_ddp_fric_b
            "Derivative at lower end of regularization domain";

          SI.Pressure dp_grav_a = g_times_height_ab*d_a
            "Static head if mass flows in design direction (a to b)";
          SI.Pressure dp_grav_b = g_times_height_ab*d_b
            "Static head if mass flows against design direction (b to a)";

          // Properly define zero mass flow conditions
          SI.MassFlowRate m_flow_zero = 0;
          SI.Pressure dp_zero = (dp_grav_a + dp_grav_b)/2;
          Real dm_flow_ddp_fric_zero;

        algorithm
          dp_a := max(dp_grav_a, dp_grav_b)+dp_small;
          dp_b := min(dp_grav_a, dp_grav_b)-dp_small;

          if dp>=dp_a then
            // Positive flow outside regularization
            m_flow := Internal.m_flow_of_dp_fric(dp-dp_grav_a, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta);
          elseif dp<=dp_b then
            // Negative flow outside regularization
            m_flow := Internal.m_flow_of_dp_fric(dp-dp_grav_b, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta);
          else
            // Regularization parameters
            (m_flow_a, dm_flow_ddp_fric_a) := Internal.m_flow_of_dp_fric(dp_a-dp_grav_a, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta);
            (m_flow_b, dm_flow_ddp_fric_b) := Internal.m_flow_of_dp_fric(dp_b-dp_grav_b, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta);
            // Include a properly defined zero mass flow point
            // Obtain a suitable slope from the linear section slope c (value of m_flow is overwritten later)
            (m_flow, dm_flow_ddp_fric_zero) := Utilities.regFun3(dp_zero, dp_b, dp_a, m_flow_b, m_flow_a, dm_flow_ddp_fric_b, dm_flow_ddp_fric_a);
            // Do regularization
            if dp>dp_zero then
              m_flow := Utilities.regFun3(dp, dp_zero, dp_a, m_flow_zero, m_flow_a, dm_flow_ddp_fric_zero, dm_flow_ddp_fric_a);
            else
              m_flow := Utilities.regFun3(dp, dp_b, dp_zero, m_flow_b, m_flow_zero, dm_flow_ddp_fric_b, dm_flow_ddp_fric_zero);
            end if;
          end if;
        end massFlowRate_dp_staticHead;

        redeclare function extends pressureLoss_m_flow_staticHead
          "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction and static head"

          annotation (smoothOrder=1);

        protected
          Real Delta = roughness/diameter "Relative roughness";
          SI.ReynoldsNumber Re1 = 745*exp(if Delta <= 0.0065 then 1 else 0.0065/Delta)
            "Boundary between laminar regime and transition";
          constant SI.ReynoldsNumber Re2 = 4000
            "Boundary between transition and turbulent regime";

          SI.MassFlowRate m_flow_a
            "Upper end of regularization domain of the dp(m_flow) relation";
          SI.MassFlowRate m_flow_b
            "Lower end of regularization domain of the dp(m_flow) relation";

          SI.Pressure dp_a "Value at upper end of regularization domain";
          SI.Pressure dp_b "Value at lower end of regularization domain";

          SI.Pressure dp_grav_a = g_times_height_ab*d_a
            "Static head if mass flows in design direction (a to b)";
          SI.Pressure dp_grav_b = g_times_height_ab*d_b
            "Static head if mass flows against design direction (b to a)";

          Real ddp_dm_flow_a
            "Derivative of pressure drop with mass flow rate at m_flow_a";
          Real ddp_dm_flow_b
            "Derivative of pressure drop with mass flow rate at m_flow_b";

          // Properly define zero mass flow conditions
          SI.MassFlowRate m_flow_zero = 0;
          SI.Pressure dp_zero = (dp_grav_a + dp_grav_b)/2;
          Real ddp_dm_flow_zero;

        algorithm
          m_flow_a := if dp_grav_a<dp_grav_b then 
            Internal.m_flow_of_dp_fric(dp_grav_b - dp_grav_a, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta)+m_flow_small else 
            m_flow_small;
          m_flow_b := if dp_grav_a<dp_grav_b then 
            Internal.m_flow_of_dp_fric(dp_grav_a - dp_grav_b, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta)-m_flow_small else 
            -m_flow_small;

          if m_flow>=m_flow_a then
            // Positive flow outside regularization
            dp := Internal.dp_fric_of_m_flow(m_flow, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta) + dp_grav_a;
          elseif m_flow<=m_flow_b then
            // Negative flow outside regularization
            dp := Internal.dp_fric_of_m_flow(m_flow, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta) + dp_grav_b;
          else
            // Regularization parameters
            (dp_a, ddp_dm_flow_a) := Internal.dp_fric_of_m_flow(m_flow_a, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta);
            dp_a := dp_a + dp_grav_a "Adding dp_grav to dp_fric to get dp";
            (dp_b, ddp_dm_flow_b) := Internal.dp_fric_of_m_flow(m_flow_b, d_a, d_b, eta_a, eta_b, length, diameter, Re1, Re2, Delta);
            dp_b := dp_b + dp_grav_b "Adding dp_grav to dp_fric to get dp";
            // Include a properly defined zero mass flow point
            // Obtain a suitable slope from the linear section slope c (value of dp is overwritten later)
            (dp, ddp_dm_flow_zero) := Utilities.regFun3(m_flow_zero, m_flow_b, m_flow_a, dp_b, dp_a, ddp_dm_flow_b, ddp_dm_flow_a);
            // Do regularization
            if m_flow>m_flow_zero then
              dp := Utilities.regFun3(m_flow, m_flow_zero, m_flow_a, dp_zero, dp_a, ddp_dm_flow_zero, ddp_dm_flow_a);
            else
              dp := Utilities.regFun3(m_flow, m_flow_b, m_flow_zero, dp_b, dp_zero, ddp_dm_flow_b, ddp_dm_flow_zero);
            end if;
          end if;
        end pressureLoss_m_flow_staticHead;

      package Internal
          "Functions to calculate mass flow rate from friction pressure drop and vice versa"
        function m_flow_of_dp_fric
            "Calculate mass flow rate as function of pressure drop due to friction"

          input SI.Pressure dp_fric
              "Pressure loss due to friction (dp = port_a.p - port_b.p)";
          input SI.Density d_a "Density at port_a";
          input SI.Density d_b "Density at port_b";
          input SI.DynamicViscosity eta_a
              "Dynamic viscosity at port_a (dummy if use_eta = false)";
          input SI.DynamicViscosity eta_b
              "Dynamic viscosity at port_b (dummy if use_eta = false)";
          input SI.Length length "Length of pipe";
          input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
          input SI.ReynoldsNumber Re1
              "Boundary between laminar regime and transition";
          input SI.ReynoldsNumber Re2
              "Boundary between transition and turbulent regime";
          input Real Delta "Relative roughness";
          output SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
          output Real dm_flow_ddp_fric
              "Derivative of mass flow rate with dp_fric";
          annotation(smoothOrder=1);

          protected
          function interpolateInRegion2_withDerivative
              "Interpolation in log-log space using a cubic Hermite polynomial, where x=log10(lambda2), y=log10(Re)"

            input Real lambda2 "Known independent variable";
            input SI.ReynoldsNumber Re1
                "Boundary between laminar regime and transition";
            input SI.ReynoldsNumber Re2
                "Boundary between transition and turbulent regime";
            input Real Delta "Relative roughness";
            input SI.Pressure dp_fric
                "Pressure loss due to friction (dp = port_a.p - port_b.p)";
            output SI.ReynoldsNumber Re "Unknown return variable";
            output Real dRe_ddp "Derivative of return value";
            annotation (smoothOrder=1);
            // point lg(lambda2(Re1)) with derivative at lg(Re1)
            protected
            Real x1=log10(64*Re1);
            Real y1=log10(Re1);
            Real y1d=1;

            // Point lg(lambda2(Re2)) with derivative at lg(Re2)
            Real aux2=Delta/3.7 + 5.74/Re2^0.9;
            Real aux3=log10(aux2);
            Real L2=0.25*(Re2/aux3)^2;
            Real aux4=2.51/sqrt(L2) + 0.27*Delta;
            Real aux5=-2*sqrt(L2)*log10(aux4);
            Real x2=log10(L2);
            Real y2=log10(aux5);
            Real y2d=0.5 + (2.51/ln(10))/(aux5*aux4);

            // Point of interest in transformed space
            Real x=log10(lambda2);
            Real y;
            Real dy_dx "Derivative in transformed space";
          algorithm
            // Interpolation
            (y, dy_dx) := Utilities.cubicHermite_withDerivative(x, x1, x2, y1, y2, y1d, y2d);

            // Return value
            Re := 10^y;

            // Derivative of return value
            dRe_ddp := Re/abs(dp_fric)*dy_dx;
          end interpolateInRegion2_withDerivative;

          SI.DynamicViscosity eta "Upstream viscosity";
          SI.Density d "Upstream density";
          Real lambda2 "Modified friction coefficient (= lambda*Re^2)";
          SI.ReynoldsNumber Re "Reynolds number";
          Real dRe_ddp "dRe/ddp";
          Real aux1;
          Real aux2;

        algorithm
          // Determine upstream density and upstream viscosity
          if dp_fric >= 0 then
            d := d_a;
            eta := eta_a;
          else
            d := d_b;
            eta := eta_b;
          end if;

          // Positive mass flow rate
          lambda2 := abs(dp_fric)*2*diameter^3*d/(length*eta*eta)
              "Known as lambda2=f(dp)";

          aux1:=(2*diameter^3*d)/(length*eta^2);

          // Determine Re and dRe/ddp under the assumption of laminar flow
          Re := lambda2/64 "Hagen-Poiseuille";
          dRe_ddp := aux1/64 "Hagen-Poiseuille";

          // Modify Re, if turbulent flow
          if Re > Re1 then
            Re :=-2*sqrt(lambda2)*log10(2.51/sqrt(lambda2) + 0.27*Delta)
                "Colebrook-White";
            aux2 := sqrt(aux1*abs(dp_fric));
            dRe_ddp := 1/ln(10)*(-2*ln(2.51/aux2+0.27*Delta)*aux1/(2*aux2)+2*2.51/(2*abs(dp_fric)*(2.51/aux2+0.27*Delta)));
            if Re < Re2 then
              (Re, dRe_ddp) := interpolateInRegion2_withDerivative(lambda2, Re1, Re2, Delta, dp_fric);
            end if;
          end if;

          // Determine mass flow rate
          m_flow := (pi*diameter/4)*eta*(if dp_fric >= 0 then Re else -Re);
          // Determine derivative of mass flow rate with dp_fric
          dm_flow_ddp_fric := (pi*diameter*eta)/4*dRe_ddp;
        end m_flow_of_dp_fric;

        function dp_fric_of_m_flow
            "Calculate pressure drop due to friction as function of mass flow rate"

          input SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
          input SI.Density d_a "Density at port_a";
          input SI.Density d_b "Density at port_b";
          input SI.DynamicViscosity eta_a
              "Dynamic viscosity at port_a (dummy if use_eta = false)";
          input SI.DynamicViscosity eta_b
              "Dynamic viscosity at port_b (dummy if use_eta = false)";
          input SI.Length length "Length of pipe";
          input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
          input SI.ReynoldsNumber Re1
              "Boundary between laminar regime and transition";
          input SI.ReynoldsNumber Re2
              "Boundary between transition and turbulent regime";
          input Real Delta "Relative roughness";
          output SI.Pressure dp_fric
              "Pressure loss due to friction (dp_fric = port_a.p - port_b.p - dp_grav)";
          output Real ddp_fric_dm_flow
              "Derivative of pressure drop with mass flow rate";
          annotation(smoothOrder=1);

          protected
          function interpolateInRegion2
              "Interpolation in log-log space using a cubic Hermite polynomial, where x=log10(Re), y=log10(lambda2)"

            input SI.ReynoldsNumber Re "Known independent variable";
            input SI.ReynoldsNumber Re1
                "Boundary between laminar regime and transition";
            input SI.ReynoldsNumber Re2
                "Boundary between transition and turbulent regime";
            input Real Delta "Relative roughness";
            input SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
            output Real lambda2 "Unknown return value";
            output Real dlambda2_dm_flow "Derivative of return value";
            annotation(smoothOrder=1);
            // point lg(lambda2(Re1)) with derivative at lg(Re1)
            protected
            Real x1 = log10(Re1);
            Real y1 = log10(64*Re1);
            Real y1d = 1;

            // Point lg(lambda2(Re2)) with derivative at lg(Re2)
            Real aux2 = Delta/3.7 + 5.74/Re2^0.9;
            Real aux3 = log10(aux2);
            Real L2 = 0.25*(Re2/aux3)^2;
            Real x2 = log10(Re2);
            Real y2 = log10(L2);
            Real y2d = 2+(2*5.74*0.9)/(ln(aux2)*Re2^0.9*aux2);

            // Point of interest in transformed space
            Real x=log10(Re);
            Real y;
            Real dy_dx "Derivative in transformed space";
          algorithm
            // Interpolation
            (y, dy_dx) := Utilities.cubicHermite_withDerivative(x, x1, x2, y1, y2, y1d, y2d);

            // Return value
            lambda2 := 10^y;

            // Derivative of return value
            dlambda2_dm_flow := lambda2/abs(m_flow)*dy_dx;
          end interpolateInRegion2;

          SI.DynamicViscosity eta "Upstream viscosity";
          SI.Density d "Upstream density";
          SI.ReynoldsNumber Re "Reynolds number";
          Real lambda2 "Modified friction coefficient (= lambda*Re^2)";
          Real dlambda2_dm_flow "dlambda2/dm_flow";
          Real aux1;
          Real aux2;

        algorithm
          // Determine upstream density and upstream viscosity
          if m_flow >= 0 then
            d := d_a;
            eta := eta_a;
          else
            d := d_b;
            eta := eta_b;
          end if;

          // Determine Reynolds number
          Re :=(4/pi)*abs(m_flow)/(diameter*eta);

          aux1 := 4/(pi*diameter*eta);

          // Use correlation for lambda2 depending on actual conditions
          if Re <= Re1 then
            lambda2 := 64*Re "Hagen-Poiseuille";
            dlambda2_dm_flow := 64*aux1 "Hagen-Poiseuille";
          elseif Re >= Re2 then
            lambda2 := 0.25*(Re/log10(Delta/3.7 + 5.74/Re^0.9))^2 "Swamee-Jain";
            aux2 := Delta/3.7+5.74/((aux1*abs(m_flow))^0.9);
            dlambda2_dm_flow := 0.5*aux1*Re*ln(10)^2*(1/(ln(aux2)^2)+(5.74*0.9)/(ln(aux2)^3*Re^0.9*aux2))
                "Swamee-Jain";
          else
            (lambda2, dlambda2_dm_flow) := interpolateInRegion2(Re, Re1, Re2, Delta, m_flow);
          end if;

          // Compute pressure drop from lambda2
          dp_fric :=length*eta*eta/(2*d*diameter*diameter*diameter)*
               (if m_flow >= 0 then lambda2 else -lambda2);

          // Compute derivative from dlambda2/dm_flow
          ddp_fric_dm_flow := (length*eta^2)/(2*diameter^3*d)*dlambda2_dm_flow;
        end dp_fric_of_m_flow;
      end Internal;
      end Detailed;

      model TestWallFrictionAndGravity
        "Pressure loss in pipe due to wall friction and gravity (only for test purposes; if needed use Pipes.StaticPipe instead)"
        extends Modelica_Fluid.Fittings.BaseClasses.PartialTwoPortTransport;

        replaceable package WallFriction = 
          Modelica_Fluid.Pipes.BaseClasses.WallFriction.QuadraticTurbulent 
          constrainedby
          Modelica_Fluid.Pipes.BaseClasses.WallFriction.PartialWallFriction
          "Characteristic of wall friction"  annotation(choicesAllMatching=true);

        parameter SI.Length length "Length of pipe";
        parameter SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
        parameter SI.Length height_ab = 0.0 "Height(port_b) - Height(port_a)" 
                                                                           annotation(Evaluate=true);
        parameter SI.Length roughness(min=0) = 2.5e-5
          "Absolute roughness of pipe (default = smooth steel pipe)" 
            annotation(Dialog(enable=WallFriction.use_roughness));

        parameter Boolean use_nominal = false
          "= true, if eta_nominal and d_nominal are used, otherwise computed from medium"
                                                                                                        annotation(Evaluate=true);
        parameter SI.DynamicViscosity eta_nominal = Medium.dynamicViscosity(
                                                       Medium.setState_pTX(
                                                           Medium.p_default, Medium.T_default, Medium.X_default))
          "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
                                                                                  annotation(Dialog(enable=use_nominal));
        parameter SI.Density d_nominal = Medium.density_pTX(Medium.p_default, Medium.T_default, Medium.X_default)
          "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
                                                                   annotation(Dialog(enable=use_nominal));

        parameter Boolean show_Re = false
          "= true, if Reynolds number is included for plotting" 
           annotation (Evaluate=true, Dialog(tab="Advanced"));
        parameter Boolean from_dp=true
          " = true, use m_flow = f(dp), otherwise dp = f(m_flow)" 
          annotation (Evaluate=true, Dialog(tab="Advanced"));
        parameter SI.AbsolutePressure dp_small = 1
          "Within regularization if |dp| < dp_small (may be wider for large discontinuities in static head)"
          annotation(Dialog(tab="Advanced", enable=from_dp and WallFriction.use_dp_small));
        parameter SI.MassFlowRate m_flow_small = reg_m_flow_small
          "Within regularizatio if |m_flow| < m_flow_small (may be wider for large discontinuities in static head)"
          annotation(Dialog(tab="Advanced", enable=not from_dp and WallFriction.use_m_flow_small));
        SI.ReynoldsNumber Re = Utilities.ReynoldsNumber_m_flow(m_flow, noEvent(if m_flow>0 then eta_a else eta_b), diameter) if show_Re
          "Reynolds number of pipe";

        outer Modelica_Fluid.System system "System properties";

      protected
        SI.DynamicViscosity eta_a = if not WallFriction.use_eta then 1.e-10 else 
                                    (if use_nominal then eta_nominal else Medium.dynamicViscosity(port_a_state_inflow));
        SI.DynamicViscosity eta_b = if not WallFriction.use_eta then 1.e-10 else 
                                    (if use_nominal then eta_nominal else Medium.dynamicViscosity(port_b_state_inflow));
        SI.Density d_a = if use_nominal then d_nominal else port_a_d_inflow;
        SI.Density d_b = if use_nominal then d_nominal else port_b_d_inflow;

        Real g_times_height_ab(final unit="m2/s2") = system.g*height_ab
          "Gravitiy times height_ab = dp_grav/d";

        // Currently not in use (means to widen the regularization domain in case of large difference in static head)
        final parameter Boolean use_x_small_staticHead = false
          "Use dp_/m_flow_small_staticHead only if static head actually exists"
                                                                                annotation(Evaluate=true);
                                                               /*abs(height_ab)>0*/
        SI.AbsolutePressure dp_small_staticHead = noEvent(max(dp_small, 0.015*abs(g_times_height_ab*(d_a-d_b))))
          "Heuristic for large discontinuities in static head";
        SI.MassFlowRate m_flow_small_staticHead = noEvent(max(m_flow_small, (-5.55e-7*(d_a+d_b)/2+5.5e-4)*abs(g_times_height_ab*(d_a-d_b))))
          "Heuristic for large discontinuities in static head";

      equation
        if from_dp and not WallFriction.dp_is_zero then
          m_flow = WallFriction.massFlowRate_dp_staticHead(dp, d_a, d_b, eta_a, eta_b, length, diameter,
            g_times_height_ab, roughness, if use_x_small_staticHead then dp_small_staticHead else dp_small);
        else
          dp = WallFriction.pressureLoss_m_flow_staticHead(m_flow, d_a, d_b, eta_a, eta_b, length, diameter,
            g_times_height_ab, roughness, if use_x_small_staticHead then m_flow_small_staticHead else m_flow_small);
        end if;

          annotation (defaultComponentName="pipeFriction",Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={1,1}), graphics={
              Rectangle(
                extent={{-100,60},{100,-60}},
                lineColor={0,0,0},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={192,192,192}),
              Rectangle(
                extent={{-100,44},{100,-45}},
                lineColor={0,0,0},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={0,127,255}),
              Text(
                extent={{-150,80},{150,120}},
                lineColor={0,0,255},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={0,127,255},
                textString="%name")}),           Documentation(info="<html>
<p>
This model describes pressure losses due to <b>wall friction</b> in a pipe
and due to gravity.
It is assumed that no mass or energy is stored in the pipe. 
Correlations of different complexity and validity can be
seleted via the replaceable package <b>WallFriction</b> (see parameter menu below).
The details of the pipe wall friction model are described in the
<a href=\"Modelica://Modelica_Fluid.UsersGuide.ComponentDefinition.WallFriction\">UsersGuide</a>.
Basically, different variants of the equation
</p>
 
<pre>
   dp = &lambda;(Re,<font face=\"Symbol\">D</font>)*(L/D)*&rho;*v*|v|/2
</pre>
 
<p>
are used, where the friction loss factor &lambda; is shown
in the next figure:
</p>
 
<img src=\"../Images/Components/PipeFriction1.png\">
 
<p>
By default, the correlations are computed with media data
at the actual time instant.
In order to reduce non-linear equation systems, parameter
<b>use_nominal</b> provides the option
to compute the correlations with constant media values
at the desired operating point. This might speed-up the
simulation and/or might give a more robust simulation.
</p>
</html>"),Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
              grid={1,1}), graphics={
              Rectangle(
                extent={{-100,64},{100,-64}},
                lineColor={0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Backward),
              Rectangle(
                extent={{-100,50},{100,-49}},
                lineColor={0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Line(
                points={{-60,-49},{-60,50}},
                color={0,0,255},
                arrow={Arrow.Filled,Arrow.Filled}),
              Text(
                extent={{-50,16},{6,-10}},
                lineColor={0,0,255},
                fillColor={0,0,255},
                fillPattern=FillPattern.Solid,
                textString="diameter"),
              Line(
                points={{-100,74},{100,74}},
                color={0,0,255},
                arrow={Arrow.Filled,Arrow.Filled}),
              Text(
                extent={{-34,92},{34,74}},
                lineColor={0,0,255},
                fillColor={0,0,255},
                fillPattern=FillPattern.Solid,
                textString="length")}));
      end TestWallFrictionAndGravity;
    end WallFriction;
  end BaseClasses;
  annotation (Documentation(info="<html>
 
</html>"));

end Pipes;
