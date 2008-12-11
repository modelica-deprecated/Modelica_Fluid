within Modelica_Fluid;
package Pipes "Lumped, distributed and thermal pipe components"
    extends Modelica_Fluid.Icons.VariantLibrary;

  model StaticPipe
    "Basic pipe flow model without storage of momentum, mass or energy"

    // extending PartialPipe
    extends Modelica_Fluid.Pipes.BaseClasses.PartialPipe(
          redeclare model HeatTransfer=BaseClasses.HeatTransfer.PipeHT_none(nPipes=nPipes));

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

    PressureDrop pressureDrop(
            redeclare final package Medium = Medium,
            final n=1,
            state={Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)),
                   Medium.setState_phX(port_b.p, inStream(port_b.h_outflow), inStream(port_b.Xi_outflow))},
            final allowFlowReversal=allowFlowReversal,
            final dynamicsType=Modelica_Fluid.Types.Dynamics.SteadyState,
            final initType=Modelica_Fluid.Types.Init.SteadyState,
            final p_a_start=p_a_start,
            final p_b_start=p_b_start,
            final m_flow_start=m_flow_start,
            final nPipes=nPipes,
            final roughness=roughness,
            final diameter=4*crossArea/perimeter,
            final length=length,
            final height_ab=height_ab,
            final g=system.g) "Pressure drop model" 
       annotation (Placement(transformation(extent={{-38,-18},{38,18}},rotation=0)));
  equation
    port_a.m_flow = pressureDrop.m_flow[1]*nPipes;
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

    // extending PartialPipe
    extends Modelica_Fluid.Pipes.BaseClasses.PartialPipe;

    // Assumptions
    parameter Modelica_Fluid.Types.Dynamics dynamicsType=system.dynamicsType
      "Dynamics option" 
      annotation(Evaluate=true, Dialog(tab = "Assumptions"));

    // Initialization
    parameter Types.Init initType=system.initType "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
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
      final nPipes=nPipes,
      final diameter=4*crossArea/perimeter,
      final area=perimeter*length,
      final crossArea=crossArea,
      final length=length,
      final use_fluidHeatPort=true,
      state={volume.medium.state},
      m_flow = {0.5*(port_a.m_flow - port_b.m_flow)}) "Heat transfer model" 
        annotation (Placement(transformation(extent={{-11,14},{11,36}},rotation=0)));

    StaticPipe staticPipe1(
      redeclare package Medium = Medium,
      allowFlowReversal=allowFlowReversal,
      nPipes=nPipes,
      length=length/2,
      roughness=roughness,
      diameter=diameter,
      perimeter=perimeter,
      crossArea=crossArea,
      height_ab=height_ab/2,
      m_flow_start=m_flow_start,
      redeclare final model PressureDrop = PressureDrop) 
      annotation (Placement(transformation(extent={{-60,-40},{-40,-20}},
            rotation=0)));
    Modelica_Fluid.Vessels.Volume volume(
      redeclare package Medium = Medium,
      initType=initType,
      p_start=(p_a_start+p_b_start)/2,
      use_T_start=use_T_start,
      T_start=T_start,
      h_start=h_start,
      X_start=X_start,
      C_start=C_start,
      dynamicsType=dynamicsType,
      V=V,
      nPorts=2,
      portDiameters={0,0},
      neglectPortDiameters=true) 
      annotation (Placement(transformation(extent={{-10,-20},{10,0}},  rotation=
             0)));
    StaticPipe staticPipe2(
      redeclare package Medium = Medium,
      allowFlowReversal=allowFlowReversal,
      nPipes=nPipes,
      length=length/2,
      roughness=roughness,
      diameter=diameter,
      perimeter=perimeter,
      crossArea=crossArea,
      height_ab=height_ab/2,
      m_flow_start=m_flow_start,
      redeclare final model PressureDrop = PressureDrop)   annotation (Placement(transformation(extent={{40,-40},
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

    // extending PartialPipe
    extends Modelica_Fluid.Pipes.BaseClasses.PartialPipe(
      final port_a_exposesState = (modelStructure == ModelStructure.av_b) or (modelStructure == ModelStructure.av_vb),
      final port_b_exposesState = (modelStructure == ModelStructure.a_vb) or (modelStructure == ModelStructure.av_vb));

    // distributed volume model
    extends Modelica_Fluid.Vessels.BaseClasses.PartialDistributedVolume(
      final n = nNodes,
      Qs_flow = heatTransfer.Q_flow);

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
      "number of flow models in pressureDrop";
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
      "state vector for pressureDrop model";

    // Pressure drop model
    PressureDrop pressureDrop(
            redeclare final package Medium = Medium,
            final n=nFlows,
            state=flowState,
            final allowFlowReversal=allowFlowReversal,
            final dynamicsType=dynamicsType,
            final initType=initType,
            final p_a_start=p_a_start,
            final p_b_start=p_b_start,
            final m_flow_start=m_flow_start,
            final nPipes=nPipes,
            final roughness=roughness,
            final diameter=4*crossArea/perimeter,
            final length=length,
            final height_ab=height_ab,
            final g=system.g) "Pressure drop model" 
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
      final nPipes=nPipes,
      final diameter=4*crossArea/perimeter,
      final area=perimeter*length,
      final crossArea=crossArea,
      final length=length,
      state=medium.state,
      m_flow = 0.5*(m_flow[1:n]+m_flow[2:n+1])) "Heat transfer model" 
        annotation (Placement(transformation(extent={{-20,-5},{20,35}},  rotation=0)));

  equation
    assert(nNodes > 1 or modelStructure <> ModelStructure.av_vb,
       "nNodes needs to be at least 2 for modelStructure av_vb, as flow model disappears otherwise!");

    // Source/sink terms for mass and energy balances
    fluidVolume=fill(V/n, n);
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
        m_flow[1] = pressureDrop.m_flow[1];
        flowState[1] = state_a;
        flowState[2] = medium[iLumped].state;
        flowState[3] = state_b;
        m_flow[n+1] = pressureDrop.m_flow[2];
      elseif modelStructure == ModelStructure.av_b then
        port_a.p = medium[1].p;
        flowState[1] = medium[iLumped].state;
        flowState[2] = state_b;
        m_flow[n+1] = pressureDrop.m_flow[1];
      elseif modelStructure == ModelStructure.a_vb then
        m_flow[1] = pressureDrop.m_flow[1];
        flowState[1] = state_a;
        flowState[2] = medium[iLumped].state;
        port_b.p = medium[n].p;
      else // av_vb
        port_a.p = medium[1].p;
        flowState[1] = medium[1].state;
        m_flow[iLumped] = pressureDrop.m_flow[1];
        flowState[2] = medium[n].state;
        port_b.p = medium[n].p;
      end if;
    else
      if modelStructure == ModelStructure.a_v_b then
        flowState[1] = state_a;
        flowState[2:n+1] = medium[1:n].state;
        flowState[n+2] = state_b;
        //m_flow = pressureDrop.m_flow;
        for i in 1:n+1 loop
          m_flow[i] = pressureDrop.m_flow[i];
        end for;
      elseif modelStructure == ModelStructure.av_b then
        flowState[1:n] = medium[1:n].state;
        flowState[n+1] = state_b;
        //m_flow[2:n+1] = pressureDrop.m_flow;
        for i in 2:n+1 loop
          m_flow[i] = pressureDrop.m_flow[i-1];
        end for;
        port_a.p = medium[1].p;
      elseif modelStructure == ModelStructure.a_vb then
        flowState[1] = state_a;
        flowState[2:n+1] = medium[1:n].state;
        //m_flow[1:n] = pressureDrop.m_flow;
        for i in 1:n loop
          m_flow[i] = pressureDrop.m_flow[i];
        end for;
        port_b.p = medium[n].p;
      else // av_vb
        flowState[1:n] = medium[1:n].state;
        //m_flow[2:n] = pressureDrop.m_flow[1:n-1];
        for i in 2:n loop
          m_flow[i] = pressureDrop.m_flow[i-1];
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
<p>Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PartialPipe\">PartialPipe</a>. 
The total volume is determined by geometry parameters. It is split into nNodes pipe segments of equal size along the flow path. 
The default value is nNodes=2.
<p><b>Mass and Energy balance</b></p>
One mass and one energy balance if formulated for each pipe segment. 
The mass and energy balances are inherited from <a href=\"Modelica:Modelica_Fluid.Vessels.BaseClasses.PartialDistributedVolume\">PartialDistributedVolume</a>. 
The additional component <b><tt>HeatTransfer</tt></b> specifies the source term <tt>Qs_flow</tt> in the energy balance. 
The default component uses a constant coefficient for the heat transfer between the bulk flow and the segment boundaries exposed through the <tt>heatPorts</tt>. 
The <tt>HeatTransfer</tt> model is replaceable and can be exchanged with any model extended from 
<a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer\">BaseClasses.PartialPipeHeatTransfer</a>.</p>
<p><b>Momentum balance</b></p>
The momentum balance is determined by the <b><tt>PressureDrop</tt></b> component, which can be replaced with any model extended from 
<a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PressureDrop.PartialPipePressureDrop\">BaseClasses.PartialPipePressureDrop</a>.
The default setting is steady-state <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PressureDrop.DetailedFlow\">DetailedFlow</a>.
The momentum balances are formed across the segment boundaries along the flow path according to the staggered grid approach. 
The default symmetric model is characterized by one momentum balance inside the pipe with nNodes=2 fluid segments.
An alternative symmetric variation with nNodes+1 momentum balances, one at each port, as well as  
non-symmetric variations can be obtained by chosing a different value for the parameter <tt><b>modelStructure</b></tt>. 
The options include:
<ul>
<li><tt>av_vb</tt>: nNodes-1 momentum balances between nNodes pipe segments, potential pressure states at both ports.
<li><tt>a_v_b</tt>: Alternative symmetric setting with nNodes+1 momentum balances across nNodes pipe segments, one momentum balance at each port. 
Connecting two pipes therefore results in algebraic pressures at the ports. 
The specification of good start values for the port pressures is essential in order to solve large systems.</li>
<li><tt>av_b</tt>: nNodes momentum balances, one between nth volume and <tt>port_b</tt>, potential pressure state at <tt>port_a</tt></li>
<li><tt>a_vb</tt>: nNodes momentum balance, one between first volume and <tt>port_a</tt>, potential pressure state at <tt>port_b</tt></li>
</ul></p>
 
<p>The PressureDrop contains
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
Consider using a junction or an adapter comonent, like <a href=\"Modelica:Modelica_Fluid.PressureLosses.SuddenExpansion\">SuddenExpansion</a>
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

    partial model PartialPipe "Base class for pipe models"
      extends Modelica_Fluid.Interfaces.PartialTwoPort;

      // Geometry

      // Note: define nPipes as Real to support inverse calculations
      parameter Real nPipes(min=1)=1 "Number of identical parallel pipes" 
        annotation(Dialog(group="Geometry"));
      parameter SI.Length length "Length"   annotation(Dialog(tab="General", group="Geometry"));
      parameter Boolean isCircular=true
        "= true if cross sectional area is circular" 
        annotation (Evaluate, Dialog(tab="General", group="Geometry"));
      parameter SI.Diameter diameter "Diameter of circular pipe"      annotation(Dialog(group="Geometry", enable=isCircular));
      parameter SI.Length perimeter=Modelica.Constants.pi*diameter
        "Inner perimeter"                                                                                       annotation(Dialog(tab="General", group="Geometry", enable=not isCircular));
      parameter SI.Area crossArea=Modelica.Constants.pi*diameter*diameter/4
        "Inner cross section area"            annotation(Dialog(tab="General", group="Geometry", enable=not isCircular));
      final parameter SI.Volume V=crossArea*length*nPipes "volume size";
      parameter SI.Length roughness(min=0)=2.5e-5
        "Average height of surface asperities (default = smooth steel pipe)" 
          annotation(Dialog(group="Geometry",enable=WallFriction.use_roughness));

      // Static head
      parameter SI.Length height_ab=0.0 "Height(port_b) - Height(port_a)" 
          annotation(Dialog(group="Static head"), Evaluate=true);

      // Pressure drop
      replaceable model PressureDrop = 
        Modelica_Fluid.Pipes.BaseClasses.PressureDrop.DetailedFlow 
        constrainedby BaseClasses.PressureDrop.PartialPipePressureDrop
        "Characteristics of wall friction and gravity" 
          annotation(Dialog(group="Pressure drop"), choicesAllMatching=true);

      // Heat transfer
      replaceable model HeatTransfer = 
          Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_constAlpha 
        constrainedby
        Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer
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

    end PartialPipe;

    package PressureDrop
      "Pressure drop models for pipes, including wall friction and static head"
          partial model PartialPipePressureDrop
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
            parameter Real nPipes "number of parallel pipes" 
               annotation(Dialog(tab="Internal Interface", enable=false,group="Geometry"));
            parameter SI.Length length "Length of flow path" 
               annotation(Dialog(tab="Internal Interface", enable=false,group="Geometry"));
            parameter SI.Diameter diameter
          "Hydraulic diameter (typically 4*crossArea/perimeter)" 
               annotation(Dialog(tab="Internal Interface", enable=false,group="Geometry"));
            parameter SI.Length roughness(min=0)
          "Average height of surface asperities" 
                annotation(Dialog(tab="Internal Interface", enable=false,group="Geometry",enable=WallFriction.use_roughness));
            parameter SI.Length height_ab
          "Height(state[n+1]) - Height(state[1])" 
                annotation(Dialog(tab="Internal Interface", enable=false,group="Static head"));

            // Additional parameters
            // Note: no outer system is used for default values,
            // as a PressureDrop model is intended as sub-component of other models
            parameter SI.Acceleration g "Constant gravity acceleration" 
              annotation(Dialog(tab="Internal Interface", enable=false,group="Static head"));
            parameter Boolean allowFlowReversal
          "= true to allow flow reversal, false restricts to design direction (state[1] -> state[n+1])"
              annotation(Dialog(tab="Internal Interface", enable=false,group="Assumptions"), Evaluate=true);
            parameter Modelica_Fluid.Types.Dynamics dynamicsType
          "Dynamics option, e.g. for models with dynamic momentum balance" 
              annotation(Dialog(tab="Internal Interface", enable=false,group = "Assumptions"), Evaluate=true);
            parameter Modelica_Fluid.Types.Init initType
          "Initialization option, e.g. for dynamic momentum balance" 
              annotation(Evaluate=true, Dialog(tab="Internal Interface", enable=false,group = "Initialization"));
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
            // Note: don't use start values for p to get same behavior as with PressureLosses.WallFrictionAndGravity
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
          end PartialPipePressureDrop;

          model NominalPressureDrop "Linear pressure drop for nominal values"
            extends PartialPipePressureDrop;

            parameter SI.AbsolutePressure dp_nominal "Nominal pressure drop";
            parameter SI.MassFlowRate m_flow_nominal
          "Mass flow rate for dp_nominal";

            parameter Boolean use_d_nominal = false
          "= true to use d_nominal for static head, otherwise computed from medium"
              annotation(Dialog(group="Advanced"), Evaluate=true);
            parameter SI.Density d_nominal = Medium.density_pTX(Medium.p_default, Medium.T_default, Medium.X_default)
          "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
              annotation(Dialog(group="Advanced", enable=use_d_nominal));

            parameter Boolean smoothFlowReversal=false
          "=true for numerical regularization around zero flow" 
              annotation(Dialog(group="Advanced",enable=allowFlowReversal and not use_d_nominal));
            parameter SI.MassFlowRate m_flow_small = 0.01
          "Within regularization if |m_flow| < m_flow_small" 
              annotation(Dialog(group="Advanced",enable=smoothFlowReversal));

            SI.Density[n+1] d = if use_d_nominal then fill(d_nominal, n+1) else Medium.density(state);

          equation
            if not allowFlowReversal or use_d_nominal then
              dp = g*height_ab*d[1:n] + dp_nominal/m_flow_nominal*m_flow*nPipes;
            else
              if not smoothFlowReversal then
                // exact switching
                for i in 1:n loop
                  dp[i] = g*height_ab*(if m_flow[i] > 0 then d[i] else d[i+1]) + dp_nominal/m_flow_nominal*m_flow[i]*nPipes;
                end for;
              else
                // regularization around zero flow
                dp = g*height_ab*Utilities.regStep(m_flow, d[1:n], d[2:n+1], m_flow_small) + dp_nominal/m_flow_nominal*m_flow*nPipes;
              end if;
            end if;
          end NominalPressureDrop;

          partial model PartialWallFrictionAndGravity
        "Base class for pressure drop in pipe due to wall friction and gravity (for both flow directions)"
            extends PartialPipePressureDrop;

            replaceable package WallFriction = 
            Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.PartialWallFriction
          "wall friction model";

            parameter Boolean use_eta_nominal = false
          "= true, if eta_nominal and d_nominal are used, otherwise computed from medium"
               annotation(Evaluate=true);
            parameter SI.DynamicViscosity eta_nominal = Medium.dynamicViscosity(
                                                           Medium.setState_pTX(
                                                               Medium.p_default, Medium.T_default, Medium.X_default))
          "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
                                                                                      annotation(Dialog(enable=use_eta_nominal));
            parameter Boolean use_d_nominal = false
          "= true, if eta_nominal and d_nominal are used, otherwise computed from medium"
               annotation(Evaluate=true);
            parameter SI.Density d_nominal = Medium.density_pTX(Medium.p_default, Medium.T_default, Medium.X_default)
          "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)"    annotation(Dialog(enable=use_d_nominal));

            parameter Boolean show_Re = false
          "= true, if Reynolds number is included for plotting" 
               annotation (Evaluate=true, Dialog(group="Advanced"));
            parameter Boolean from_dp=true
          " = true, use m_flow = f(dp), otherwise dp = f(m_flow)" 
              annotation (Evaluate=true, Dialog(group="Advanced"));
            parameter SI.AbsolutePressure dp_small = 1
          "Within regularization if |dp| < dp_small (may be wider for large discontinuities in static head)"
              annotation(Dialog(group="Advanced", enable=from_dp and WallFriction.use_dp_small));
            parameter SI.MassFlowRate m_flow_small = 0.01
          "Within regularization if |m_flow| < m_flow_small (may be wider for large discontinuities in static head)"
              annotation(Dialog(group="Advanced", enable=not from_dp and WallFriction.use_m_flow_small));
            SI.ReynoldsNumber[n] Re=Modelica_Fluid.Utilities.ReynoldsNumber_m_flow(
                m_flow/nPipes,
                (eta[1:n] + eta[2:n+1])*0.5,
                diameter) if show_Re "Reynolds numbers of pipe flow";

            // internal variables
            SI.DynamicViscosity[n+1] eta = if not WallFriction.use_eta then fill(1e-10, n+1) else 
                                        (if use_eta_nominal then fill(eta_nominal, n+1) else Medium.dynamicViscosity(state));
            SI.Density[n+1] d = if use_d_nominal then fill(d_nominal, n+1) else Medium.density(state);

            Real g_times_height_ab(final unit="m2/s2") = g*height_ab
          "Gravitiy times height_ab = dp_grav/d";

            // basing static head on average density is sensible for a distributed pipe with small segments
            final parameter Boolean use_staticHead = true
          "= false to use average density for static head, independent of regularization of flow reversal"
                                                                                               annotation(Dialog(group="Advanced"), Evaluate=true);

          equation
          if use_staticHead then
            // regularization with static head
            if from_dp and not WallFriction.dp_is_zero then
              m_flow = WallFriction.massFlowRate_dp_staticHead(
                dp,
                d[1:n],
                d[2:n+1],
                eta[1:n],
                eta[2:n+1],
                length/n,
                diameter,
                g_times_height_ab/n,
                roughness,
                dp_small/n)*nPipes;
            else
              dp = WallFriction.pressureLoss_m_flow_staticHead(
                m_flow/nPipes,
                d[1:n],
                d[2:n+1],
                eta[1:n],
                eta[2:n+1],
                length/n,
                diameter,
                g_times_height_ab/n,
                roughness,
                m_flow_small/nPipes);
            end if;
          else
            if from_dp and not WallFriction.dp_is_zero then
              m_flow = WallFriction.massFlowRate_dp(
                dp - g*height_ab/n*(d[1:n] + d[2:n+1])/2,
                d[1:n],
                d[2:n+1],
                eta[1:n],
                eta[2:n+1],
                length/n,
                diameter,
                roughness,
                dp_small)*nPipes;
            else
              dp = WallFriction.pressureLoss_m_flow(
                m_flow/nPipes,
                d[1:n],
                d[2:n+1],
                eta[1:n],
                eta[2:n+1],
                length/n,
                diameter,
                roughness,
                m_flow_small/nPipes) + g*height_ab/n*(d[1:n] + d[2:n+1])/2;
            end if;
          end if;

              annotation (defaultComponentName="pipeFriction",Icon(coordinateSystem(
                  preserveAspectRatio=false,
                  extent={{-100,-100},{100,100}},
                  grid={1,1}), graphics={Text(
                extent={{-150,80},{150,120}},
                lineColor={0,0,255},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={0,127,255},
                textString="%name")}),               Documentation(info="<html>
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
          end PartialWallFrictionAndGravity;

          partial model NoWallFriction
        "No pipe wall friction ... partial as height_ab is not supported"
           extends PartialWallFrictionAndGravity(
              redeclare package WallFriction = 
              PressureLosses.BaseClasses.WallFriction.NoFriction,
              final from_dp = false);
          equation
            assert(abs(height_ab) < Modelica.Constants.small,
             "Current implementation of WallFriction.NoFriction (r1953) does not consider static head!");
          end NoWallFriction;

          model LaminarFlow
        "Pipe wall friction in the laminar regime (linear correlation)"
           extends PartialWallFrictionAndGravity(
              redeclare package WallFriction = 
              PressureLosses.BaseClasses.WallFriction.Laminar);
          end LaminarFlow;

          model QuadraticTurbulentFlow
        "Pipe wall friction in the quadratic turbulent regime (simple characteristic, eta not used)"
           extends PartialWallFrictionAndGravity(
              redeclare package WallFriction = 
              PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent);
          end QuadraticTurbulentFlow;

          model LaminarAndQuadraticTurbulentFlow
        "Pipe wall friction in the laminar and quadratic turbulent regime (simple characteristic)"
           extends PartialWallFrictionAndGravity(
              redeclare package WallFriction = 
              PressureLosses.BaseClasses.WallFriction.LaminarAndQuadraticTurbulent);
          end LaminarAndQuadraticTurbulentFlow;

          model DetailedFlow
        "Pipe wall friction in the whole regime (detailed characteristic)"
           extends PartialWallFrictionAndGravity(
              redeclare package WallFriction = 
              PressureLosses.BaseClasses.WallFriction.Detailed);
          end DetailedFlow;
    end PressureDrop;

  package HeatTransfer
    partial model PartialPipeHeatTransfer
        "base class for any pipe heat transfer correlation"

      // Parameters
      replaceable package Medium=Modelica.Media.Interfaces.PartialMedium 
        annotation(Dialog(tab="Internal Interface", enable=false));
      parameter Integer n=1 "Number of heat transfer segments" 
        annotation(Dialog(tab="Internal Interface", enable=false), Evaluate=true);
      parameter Real nPipes "Number of parallel pipes" 
        annotation(Dialog(tab="Internal Interface", enable=false));
      parameter SI.Length length "Pipe length" 
        annotation(Dialog(tab="Internal Interface", enable=false));
      parameter SI.Length diameter
          "Hydraulic diameter (typically 4*crossArea/perimeter)" 
        annotation(Dialog(tab="Internal Interface", enable=false));
      parameter SI.Area crossArea "Cross flow area" 
        annotation(Dialog(tab="Internal Interface", enable=false));
      parameter SI.Area area "Heat transfer area (typically perimeter*length)" 
        annotation(Dialog(tab="Internal Interface", enable=false));
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
    end PartialPipeHeatTransfer;

    partial model PartialPipeHT_Nu
        "Base class for pipe heat transfer correlation in terms of Nusselt numberheat transfer in a circular pipe for laminar and turbulent one-phase flow"
      extends PartialPipeHeatTransfer;
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
      Re = CharacteristicNumbers.ReynoldsNumber(m_flow/nPipes, diameter, crossArea, eta);
      Nu = CharacteristicNumbers.NusseltNumber(alpha, diameter, lambda);
      Q_flow={alpha[i]*area/n*(wallHeatPort[i].T - T[i])*nPipes for i in 1:n};
        annotation (Documentation(info="<html>
Base class for heat transfer models that are expressed in terms of the Nusselt number and which can be used in distributed pipe models.
</html>"));
    end PartialPipeHT_Nu;

    model PipeHT_none
        "PipeHT_none: No heat transfer assuming perfect isolation"
      extends PartialPipeHeatTransfer;
    equation
      Q_flow = zeros(n);
      annotation(Documentation(info="<html>
Ideal heat transfer without thermal resistance.
</html>"));
    end PipeHT_none;

    model PipeHT_ideal
        "PipeHT_ideal: Ideal heat transfer without thermal resistance"
      extends PartialPipeHeatTransfer;
    equation
      T = wallHeatPort.T;
      annotation(Documentation(info="<html>
Ideal heat transfer without thermal resistance.
</html>"));
    end PipeHT_ideal;

    model PipeHT_constAlpha
        "PipeHT_constAlpha: Constant heat transfer coefficient"
      extends PartialPipeHeatTransfer;
      parameter SI.CoefficientOfHeatTransfer alpha0=200
          "heat transfer coefficient";
      annotation(Documentation(info="<html>
Simple heat transfer correlation with constant heat transfer coefficient, used as default component in <a distributed pipe models.
</html>"));
    equation
      Q_flow = alpha0*area/n*(wallHeatPort.T - T)*nPipes;
    end PipeHT_constAlpha;
    annotation (Documentation(info="<html>
Heat transfer correlations for pipe models
</html>"));

    model PipeHT_localLamTurb
        "PipeHT_localLamTurb: Laminar and turbulent forced convection in pipes, local coefficients"
      extends PartialPipeHT_Nu;
      protected
      Real[n] Nu_turb "Nusselt number for turbulent flow";
      Real[n] Nu_lam "Nusselt number for laminar flow";
      Real Nu_1;
      Real[n] Nu_2;
      Real[n] Xi;
      parameter SI.Length dx=length/n;
    equation
      Nu_1=3.66;
      for i in 1:n loop
       Nu_turb[i]=smooth(0,(Xi[i]/8)*abs(Re[i])*Pr[i]/(1+12.7*(Xi[i]/8)^0.5*(Pr[i]^(2/3)-1))*(1+1/3*(diameter/dx/(if m_flow[i]>=0 then (i-0.5) else (n-i+0.5)))^(2/3)));
       Xi[i]=(1.8*Modelica.Math.log10(max(1e-10,Re[i]))-1.5)^(-2);
       Nu_lam[i]=(Nu_1^3+0.7^3+(Nu_2[i]-0.7)^3)^(1/3);
       Nu_2[i]=smooth(0,1.077*(abs(Re[i])*Pr[i]*diameter/dx/(if m_flow[i]>=0 then (i-0.5) else (n-i+0.5)))^(1/3));
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
    end PipeHT_localLamTurb;
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
  end BaseClasses;
  annotation (Documentation(info="<html>
 
</html>"));

end Pipes;
