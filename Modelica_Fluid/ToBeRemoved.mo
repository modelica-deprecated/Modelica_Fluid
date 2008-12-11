within Modelica_Fluid;
package ToBeRemoved "models that will disappear from the release"
  model StaticPipe_Old
    "Basic pipe flow model without storage of mass or energy"
    extends Modelica_Fluid.ToBeRemoved.PartialPipe_Old(
       redeclare model HeatTransfer = 
          Pipes.BaseClasses.HeatTransfer.PipeHT_ideal);
    Modelica_Fluid.PressureLosses.WallFrictionAndGravity wallFriction(
      redeclare package Medium = Medium,
      allowFlowReversal=allowFlowReversal,
      redeclare package WallFriction = WallFriction,
      roughness=roughness,
      eta_nominal=eta_nominal,
      d_nominal=d_nominal,
      dp_small=dp_small,
      m_flow_start=m_flow_start,
      compute_T=false,
      show_Re=show_Re,
      from_dp=from_dp,
      diameter=diameter_h,
      reg_m_flow_small=m_flow_small,
      m_flow_small=m_flow_small,
      dp_start=(p_a_start - p_b_start),
      length=length,
      height_ab=height_ab,
      use_nominal=use_eta_nominal or use_d_nominal) 
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
            rotation=0)));
  equation
    connect(port_a, wallFriction.port_a) annotation (Line(
        points={{-100,0},{-10,0}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(wallFriction.port_b, port_b) annotation (Line(
        points={{10,0},{100,0}},
        color={0,127,255},
        smooth=Smooth.None));
    annotation (defaultComponentName="pipe", Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
              -100},{100,100}}),      graphics));
  end StaticPipe_Old;

  model LumpedPipe_Old
    // Assumptions
    parameter Modelica_Fluid.Types.Dynamics dynamicsType=system.dynamicsType
      "Dynamics option" 
      annotation(Evaluate=true, Dialog(tab = "Assumptions"));

    // Initialization
    parameter Types.Init initType=system.initType "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));

    // Extend here to get right ordering in parameter box
    extends Modelica_Fluid.ToBeRemoved.PartialPipe_Old;

     //Initialization
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

    replaceable HeatTransfer heatTransfer(
      redeclare final package Medium = Medium,
      final n=1,
      final nPipes=1,
      diameter=4*crossArea/perimeter,
      area=perimeter*length,
      final crossArea=crossArea,
      final length=length,
      state={volume.medium.state},
      m_flow = {0.5*(port_a.m_flow - port_b.m_flow)},
      final use_fluidHeatPort=true) "Edit heat transfer parameters" 
      annotation (editButton=true, Placement(transformation(extent={{-20,0},{20,40}},   rotation=0)));

    Modelica_Fluid.PressureLosses.WallFrictionAndGravity wallFriction1(
      redeclare package Medium = Medium,
      allowFlowReversal=allowFlowReversal,
      redeclare package WallFriction = WallFriction,
      length=length/2,
      height_ab=height_ab/2,
      roughness=roughness,
      eta_nominal=eta_nominal,
      d_nominal=d_nominal,
      dp_small=dp_small,
      dp_start = (p_a_start - p_b_start)/2,
      m_flow_start = m_flow_start,
      compute_T=false,
      show_Re=show_Re,
      from_dp=from_dp,
      diameter=diameter_h,
      reg_m_flow_small=m_flow_small,
      m_flow_small=m_flow_small,
      use_nominal=use_eta_nominal or use_d_nominal) 
      annotation (Placement(transformation(extent={{-60,-30},{-40,-10}},
            rotation=0)));
    Modelica_Fluid.Vessels.Volume volume(
      redeclare package Medium = Medium,
      initType=initType,
      p_start=(p_a_start+p_b_start)/2,
      use_T_start=use_T_start,
      T_start=T_start,
      h_start=h_start,
      X_start=X_start,
      dynamicsType=dynamicsType,
      V=V,
      nPorts=2,
      portDiameters={0,0},
      neglectPortDiameters=true) 
      annotation (Placement(transformation(extent={{-10,-30},{10,-10}},rotation=
             0)));
    Modelica_Fluid.PressureLosses.WallFrictionAndGravity wallFriction2(
      redeclare package Medium = Medium,
      allowFlowReversal=allowFlowReversal,
      redeclare package WallFriction = WallFriction,
      length=length/2,
      height_ab=height_ab/2,
      roughness=roughness,
      eta_nominal=eta_nominal,
      d_nominal=d_nominal,
      dp_small=dp_small,
      dp_start = (p_a_start - p_b_start)/2,
      m_flow_start = m_flow_start,
      compute_T=false,
      show_Re=show_Re,
      from_dp=from_dp,
      m_flow_small=m_flow_small,
      diameter=diameter_h,
      reg_m_flow_small=m_flow_small,
      use_nominal=use_eta_nominal or use_d_nominal) 
                                    annotation (Placement(transformation(extent={{40,-30},
              {60,-10}},         rotation=0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort 
      annotation (Placement(transformation(extent={{-10,44},{10,64}}, rotation=
              0)));
  equation
    connect(wallFriction1.port_a, port_a) 
      annotation (Line(points={{-60,-20},{-80,-20},{-80,0},{-100,0}},
                                                  color={0,127,255}));
    connect(wallFriction2.port_b, port_b) 
      annotation (Line(points={{60,-20},{80,-20},{80,0},{100,0}},
                                                color={0,127,255}));
    connect(heatPort, heatTransfer.wallHeatPort[1]) annotation (Line(
        points={{0,54},{0,34}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(heatTransfer.fluidHeatPort[1], volume.heatPort) annotation (Line(
        points={{0,8},{0,-10.2}},
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
    connect(wallFriction1.port_b, volume.ports[1]) annotation (Line(
        points={{-40,-20},{-25,-20},{-25,-40},{0,-40},{0,-28}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(wallFriction2.port_a, volume.ports[2]) annotation (Line(
        points={{40,-20},{24,-20},{24,-40},{0,-40},{0,-32}},
        color={0,127,255},
        smooth=Smooth.None));
  end LumpedPipe_Old;

 model DistributedPipeLumpedPressure_Old
    "Distributed pipe model with lumped pressure state"
    import Modelica_Fluid.Types.ModelStructure;
   extends Modelica_Fluid.ToBeRemoved.PartialDistributedFlow_Old(
     Qs_flow=heatTransfer.Q_flow,
     final port_a_exposesState = (modelStructure == ModelStructure.av_b) or (modelStructure == ModelStructure.av_vb),
     final port_b_exposesState = (modelStructure == ModelStructure.a_vb) or (modelStructure == ModelStructure.av_vb));

   parameter Types.ModelStructure modelStructure=Types.ModelStructure.a_v_b
      "Determines whether flow or volume models are present at the ports"                             annotation(Evaluate=true);

   final parameter Integer nl=integer(n/2)+1
      "Number of control volume that contains single pressure state"                  annotation(Evaluate=true);

   final parameter SI.Pressure[2] dp_start = {p_a_start - p_start[nl], p_start[nl] - p_b_start};
   SI.Pressure[2] dp(start=dp_start)
      "Pressure difference across staggered grid";

   SI.ReynoldsNumber[n+1] Re=Modelica_Fluid.Utilities.ReynoldsNumber_m_flow(
       m_flow,
       (cat(1, {eta_a}, eta) + cat(1, eta, {eta_b}))*0.5,
       diameter) if                                                                                              show_Re
      "Reynolds number of pipe flow";
   HeatTransfer heatTransfer(
     redeclare final package Medium = Medium,
     final n=nNodes,
     final nPipes=1,
     diameter=4*crossArea/perimeter,
     area=perimeter*length,
     final crossArea=crossArea,
     final length=length,
     state=medium.state,
     m_flow = 0.5*(m_flow[1:n]+m_flow[2:n+1])) "Edit heat transfer parameters" 
             annotation (Placement(transformation(extent={{-20,-5},{20,35}},  rotation=0)));

   SI.Length[2] dlength "discretized length for pressure drop";
   SI.Length[2] dheight_ab "discretized height_ab for static head";

  protected
   SI.DynamicViscosity eta_a=if not WallFriction.use_eta then 1.e-10 else (if 
       use_eta_nominal then eta_nominal else (if useInnerPortProperties then Medium.dynamicViscosity(medium[1].state) else Medium.dynamicViscosity(Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)))));
   SI.DynamicViscosity eta_b=if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then eta_nominal else (if useInnerPortProperties then Medium.dynamicViscosity(medium[n].state) else Medium.dynamicViscosity(Medium.setState_phX(port_b.p, inStream(port_b.h_outflow), inStream(port_b.Xi_outflow)))));

 equation
   // Only one connection allowed to a port to avoid unwanted ideal mixing
   assert(cardinality(port_a) <= 1 or (modelStructure == ModelStructure.a_vb) or (modelStructure == ModelStructure.a_v_b),"
port_a exposing volume with selected modelStructure shall at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller. Use a Junctions.MultiPort.
");
   assert(cardinality(port_b) <= 1 or (modelStructure == ModelStructure.av_b) or (modelStructure == ModelStructure.a_v_b),"
port_b exposing volume with selected modelStructure shall at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller. Use a Junctions.MultiPort.
");

   Ws_flow=zeros(n);
   fluidVolume=ones(n)*V/n;

   //Momentum Balance, dp contains contributions from acceleration, gravitational and friction effects
   //two momentum balances, one on each side of pressure state
   dp = {port_a.p - p[nl], p[nl] - port_b.p};
   //lumped pressure
   p[1]*ones(n-1) = p[2:n];

   if modelStructure == ModelStructure.a_v_b then
     dlength[1] = ((integer(n/2) + 1)*2 - 1)/(2*n)*length;
     dlength[2] = (2*n - (integer(n/2) + 1)*2 + 1)/(2*n)*length;
   elseif modelStructure == ModelStructure.av_b then
     dlength = {0, length};
   elseif modelStructure == ModelStructure.a_vb then
     dlength = {length, 0};
   else // av_vb
     dlength = {0, 0};
   end if;
   dheight_ab = dlength/length*height_ab;

  if from_dp and not WallFriction.dp_is_zero then
    if port_a_exposesState then
      port_a.p = p[1];
    else
      m_flow[1] = WallFriction.massFlowRate_dp_staticHead(
        dp[1],
        d_a,
        d_b,
        eta_a,
        eta_b,
        dlength[1],
        diameter_h,
        dheight_ab[1]*system.g,
        roughness,
        dp_small);
    end if;
    if port_b_exposesState then
      port_b.p = p[n];
    else
      m_flow[n + 1] = WallFriction.massFlowRate_dp_staticHead(
        dp[2],
        d_a,
        d_b,
        eta_a,
        eta_b,
        dlength[2],
        diameter_h,
        dheight_ab[2]*system.g,
        roughness,
        dp_small);
     end if;
  else
    if port_a_exposesState then
      port_a.p = p[1];
    else
      dp[1] = WallFriction.pressureLoss_m_flow_staticHead(
        m_flow[1],
        d_a,
        d_b,
        eta_a,
        eta_b,
        dlength[1],
        diameter_h,
        dheight_ab[1]*system.g,
        roughness,
        m_flow_small);
    end if;
    if port_b_exposesState then
      port_b.p = p[n];
    else
      dp[2] = WallFriction.pressureLoss_m_flow_staticHead(
        m_flow[n+1],
        d_a,
        d_b,
        eta_a,
        eta_b,
        dlength[2],
        diameter_h,
        dheight_ab[2]*system.g,
        roughness,
        m_flow_small);
    end if;
  end if;

   connect(heatPorts, heatTransfer.wallHeatPort) 
     annotation (Line(points={{0,54},{0,29}}, color={191,0,0}));
   annotation (defaultComponentName="pipe",
     Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}},
          grid={1,1}), graphics={
          Ellipse(
            extent={{-72,10},{-52,-10}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-30,10},{-10,-10}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{10,10},{30,-10}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{50,10},{70,-10}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid)}),
     Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
              {100,100}},
          grid={1,1}),
             graphics),
     Documentation(info="<html>
Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow_pLumped\">PartialDistributedFlowLumpedPressure</a>. Source terms in mass and energy balances are set to zero. The total volume is a paramter. The number of momentum balances is reduced to two, one on each side of the hydraulic state, which corresponds to a constant pressure along the entire pipe with pressure drop and gravitational forces lumped at the ports.<The additional component <tt>heatTransfer</tt> specifies the source term <tt>Qs_flow</tt> in the energy balance. The default component uses a constant coefficient of heat transfer to model convective heat transfer between segment boundary (<tt>heatPorts</tt>) and the bulk flow. The <tt>heatTransfer</tt> model is replaceable and can be exchanged with any model extended from <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer\">PartialPipeHeatTransfer</a>. .
<p><b>Momentum balance</b></p>
<p>The momentum balance is always static, i.e. no dynamic momentum term is used. The momentum balances are formed across the segment boundaries (staggered grid). For this model only two momentum balances are formed on each side of a single pressure state (roughly half way along the flowpath). This assumes a constant pressure level for all medium models in the pipe. The total pressure drop (or rise) is split to be located on each end of the component. Connecting two pipes results in an algebraic pressure at the ports. Specifying a good start value for the port pressure is essential in order to solve large systems. The term <tt>dp</tt> is unspecified in this partial class. When extending from this model it may contain
<ul>
<li>pressure drop due to friction and other dissipative losses</li>
<li>gravity effects for non-horizontal pipes</li>
</ul>
Not considered are changes in pressure resulting from significant variation of flow velocity along the flow path (with the assumption of a constant cross sectional area it must result from fluid density changes, such as in two-phase flow)</li>
 
When connecting two components, e.g. two pipes, the momentum balance across the connection point reduces to</p> 
<pre>pipe1.port_b.p = pipe2.port_a.p</pre>
<p>This is only true if the flow velocity remains the same on each side of the connection. For any significant change in diameter (and if the resulting effects, such as change in kinetic energy, cannot be neglected) an adapter component should be used. This also allows for taking into account friction losses with respect to the actual geometry of the connection point.</p>
</html>",  revisions="<html>
<ul>
<li><i>27 Nov 2008</i>
    by R&uuml;diger Franke:<br>
       Major overhaul</li>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));

 end DistributedPipeLumpedPressure_Old;

  model DistributedPipe_Old "Distributed pipe model"
    import Modelica_Fluid.Types.ModelStructure;
  extends Modelica_Fluid.ToBeRemoved.PartialDistributedFlow_Old(
    Qs_flow=heatTransfer.Q_flow,
    final port_a_exposesState = (modelStructure == ModelStructure.av_b) or (modelStructure == ModelStructure.av_vb),
    final port_b_exposesState = (modelStructure == ModelStructure.a_vb) or (modelStructure == ModelStructure.av_vb));

    parameter Types.ModelStructure modelStructure=Types.ModelStructure.a_v_b
      "Determines whether flow or volume models are present at the ports"                              annotation(Evaluate=true);

    parameter Boolean use_staticHead = true
      "= true, if WallFriction.*_staticHead functions should be used"                          annotation(Dialog(tab="Advanced"),Evaluate=true);

    // pressure loss and static head
    final parameter SI.Pressure dp_start[n+1]=cat(1, {p_a_start}, p_start) - cat(1, p_start, {p_b_start});
    SI.Pressure[n+1] dp(start=dp_start)
      "pressure difference across staggered grid";

  SI.ReynoldsNumber[n+1] Re=Modelica_Fluid.Utilities.ReynoldsNumber_m_flow(
    m_flow,
    (cat(1, {eta_a}, eta) + cat(1, eta, {eta_b}))*0.5,
    diameter) if                                                                                               show_Re
      "Reynolds number of pipe flow";
  replaceable HeatTransfer heatTransfer(
    redeclare each final package Medium = Medium,
    final n=n,
    final nPipes=1,
    diameter=4*crossArea/perimeter,
    area=perimeter*length,
    final crossArea=crossArea,
    final length=length,
    state=medium.state,
    m_flow = 0.5*(m_flow[1:n]+m_flow[2:n+1])) "Edit heat transfer parameters" 
              annotation (editButton=true, Placement(transformation(extent={{-20,-5},{20,35}},  rotation=0)));

  SI.Length[n+1] dlength "discretized length for pressure drop";
  SI.Length[n+1] dheight_ab "discretized height_ab for static head";
  protected
  SI.DynamicViscosity eta_a=if not WallFriction.use_eta then 1.e-10 else (if 
      use_eta_nominal then eta_nominal else (if useInnerPortProperties then eta[1] else Medium.dynamicViscosity(Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)))));
  SI.DynamicViscosity eta_b=if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then eta_nominal else (if useInnerPortProperties then eta[n] else Medium.dynamicViscosity(Medium.setState_phX(port_b.p, inStream(port_b.h_outflow), inStream(port_b.Xi_outflow)))));
  SI.DynamicViscosity[n] eta=if not WallFriction.use_eta then fill(1.e-10, n) else (if use_eta_nominal then fill(eta_nominal, n) else 
      Medium.dynamicViscosity(medium.state));

  equation
    // Only one connection allowed to a port to avoid unwanted ideal mixing
    assert(cardinality(port_a) <= 1 or (modelStructure == ModelStructure.a_vb) or (modelStructure == ModelStructure.a_v_b),"
port_a exposing volume with selected modelStructure shall at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller. Use a Junctions.MultiPort.
");
    assert(cardinality(port_b) <= 1 or (modelStructure == ModelStructure.av_b) or (modelStructure == ModelStructure.a_v_b),"
port_b exposing volume with selected modelStructure shall at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller. Use a Junctions.MultiPort.
");

    Ws_flow=zeros(n);
    fluidVolume=ones(n)*V/n;

    //Pressure drop and gravity
    //Simplified Momentum Balance, dp contains contributions from gravitational and friction effects
    dp = cat(1, {port_a.p}, p) - cat(1, p, {port_b.p});

    if modelStructure == ModelStructure.a_v_b then
      dlength = cat(1, {0.5*length/n}, length/n*ones(n-1), {0.5*length/n});
    elseif modelStructure == ModelStructure.av_b then
      dlength = cat(1, {0}, length/n*ones(n));
    elseif modelStructure == ModelStructure.a_vb then
      dlength = cat(1, length/n*ones(n), {0});
    else // av_vb
      dlength = cat(1, {0}, length/(n-1)*ones(n-1), {0});
    end if;
    dheight_ab = dlength/length*height_ab;

    //Pressure drop and gravity
  if use_staticHead then
    if from_dp and not WallFriction.dp_is_zero then
      if port_a_exposesState then
        port_a.p = medium[1].p;
      else
        m_flow[1] = WallFriction.massFlowRate_dp_staticHead(
          dp[1],
          d_a,
          d[1],
          eta_a,
          eta[1],
          dlength[1],
          diameter_h,
          dheight_ab[1]*system.g,
          roughness,
          dp_small);
      end if;
      for i in 2:n loop
        m_flow[i] = WallFriction.massFlowRate_dp_staticHead(
          dp[i],
          d[i - 1],
          d[i],
          eta[i - 1],
          eta[i],
          dlength[i],
          diameter_h,
          dheight_ab[i]*system.g,
          roughness,
          dp_small);
      end for;
      if port_b_exposesState then
        port_b.p=medium[n].p;
      else
        m_flow[n + 1] = WallFriction.massFlowRate_dp_staticHead(
          dp[n+1],
          d[n],
          d_b,
          eta[n],
          eta_b,
          dlength[n+1],
          diameter_h,
          dheight_ab[n+1]*system.g,
          roughness,
          dp_small);
      end if;

    else

      if port_a_exposesState then
        port_a.p = medium[1].p;
      else
        dp[1] = WallFriction.pressureLoss_m_flow_staticHead(
          m_flow[1],
          d_a,
          d[1],
          eta_a,
          eta[1],
          dlength[1],
          diameter_h,
          dheight_ab[1]*system.g,
          roughness,
          m_flow_small);
      end if;
      for i in 2:n loop
        dp[i] = WallFriction.pressureLoss_m_flow_staticHead(
          m_flow[i],
          d[i - 1],
          d[i],
          eta[i - 1],
          eta[i],
          dlength[i],
          diameter_h,
          dheight_ab[i]*system.g,
          roughness,
          m_flow_small);
      end for;
      if port_b_exposesState then
        port_b.p=medium[n].p;
      else
        dp[n+1] = WallFriction.pressureLoss_m_flow_staticHead(
          m_flow[n+1],
          d[n],
          d_b,
          eta[n],
          eta_b,
          dlength[n+1],
          diameter_h,
          dheight_ab[n+1]*system.g,
          roughness,
          m_flow_small);
      end if;
    end if;

  else // no *_staticHead

    if from_dp and not WallFriction.dp_is_zero then

    if port_a_exposesState then
      port_a.p = medium[1].p;
    else
      m_flow[1] = WallFriction.massFlowRate_dp(
        dp[1] - dheight_ab[1]*system.g*d[1],
        d_a,
        d[1],
        eta_a,
        eta[1],
        dlength[1],
        diameter_h,
        roughness,
        dp_small);
    end if;
    for i in 2:n loop
      m_flow[i] = WallFriction.massFlowRate_dp(
        dp[i] - dheight_ab[i]*system.g*(d[i-1] + d[i])/2,
        d[i - 1],
        d[i],
        eta[i - 1],
        eta[i],
        dlength[i],
        diameter_h,
        roughness,
        dp_small);
    end for;
    if port_b_exposesState then
      port_b.p=medium[n].p;
    else
      m_flow[n + 1] = WallFriction.massFlowRate_dp(
        dp[n+1] - dheight_ab[n+1]*system.g*d[n],
        d[n],
        d_b,
        eta[n],
        eta_b,
        dlength[n + 1],
        diameter_h,
        roughness,
        dp_small);
    end if;

    else

    if port_a_exposesState then
      port_a.p = medium[1].p;
    else
      dp[1] = WallFriction.pressureLoss_m_flow(
        m_flow[1],
        d_a,
        d[1],
        eta_a,
        eta[1],
        dlength[1],
        diameter_h,
        roughness,
        m_flow_small) + height_ab/n*system.g*d[1]/2;
    end if;
    for i in 2:n loop
      dp[i] = WallFriction.pressureLoss_m_flow(
        m_flow[i],
        d[i - 1],
        d[i],
        eta[i - 1],
        eta[i],
        dlength[i],
        diameter_h,
        roughness,
        m_flow_small) + height_ab/n*system.g*(d[i - 1] + d[i])/2;
    end for;
    if port_b_exposesState then
      port_b.p=medium[n].p;
    else
      dp[n+1] = WallFriction.pressureLoss_m_flow(
        m_flow[n+1],
        d[n],
        d_b,
        eta[n],
        eta_b,
        dlength[n+1],
        diameter_h,
        roughness,
        m_flow_small) + height_ab/n/2*system.g*d[n];
    end if;
    end if;

  end if;

    connect(heatPorts, heatTransfer.wallHeatPort) 
      annotation (Line(points={{0,54},{0,29}}, color={191,0,0}));

    annotation (defaultComponentName="pipe",
  Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={
          Ellipse(
            extent={{-72,10},{-52,-10}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-30,10},{-10,-10}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{10,10},{30,-10}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{50,10},{70,-10}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid)}),
  Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{100,
              100}},
          grid={1,1}),
          graphics),
  Documentation(info="<html>
<p>Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow\">PartialDistributedFlow</a>. Source terms in the mass balances are set to zero. The total volume is a parameter. The additional component <tt>heatTransfer</tt> specifies the source term <tt>Qs_flow</tt> in the energy balance. The default component uses a constant coefficient of heat transfer to model convective heat transfer between segment boundary (<tt>heatPorts</tt>) and the bulk flow. The <tt>heatTransfer</tt> model is replaceable and can be exchanged with any model extended from <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer\">PartialPipeHeatTransfer</a>.</p>
<p>Pressure drop correlations (algebraic and possibly non-linear flow model) correlate the pressure in the first control volume with the pressure in port_a and the pressures of port_b and the nth control volume, respectively.</p>
<p><b>Momentum balance</b></p>
<p>The momentum balance is always static, i.e. no dynamic momentum term is used. The momentum balances are formed across the segment boundaries (staggered grid). The default symmetric model is characterized by half a momentum balance on each end of the flow model resulting in a total of n-1 full and 2 half momentum balances. Connecting two pipes therefore results in an algebraic pressure at the ports. Specifying a good start value for the port pressure is essential in order to solve large systems. Non-symmetric variations are obtained by chosing a different value for the parameter <tt><b>modelStructure</b></tt>. Options include:
<ul>
<li><tt>a_v_b</tt>: default setting with two half momentum balances</li>
<li><tt>av_b</tt>: full momentum balance between nth volume and <tt>port_b</tt>, potential pressure state at <tt>port_a</tt></li>
<li><tt>a_vb</tt>: full momentum balance between first volume and <tt>port_a</tt>, potential pressure state at <tt>port_b</tt></li>
<li><tt>av_vb</tt>: n-1 momentum balances between first and nth volume, potential pressure states at both ports. It's use should be avoided, since not the entire pipe length is taken into account.
</ul></p>
 
<p>The term <tt>dp</tt> contains
<ul>
<li>pressure drop due to friction and other dissipative losses</li>
<li>gravity effects for non-horizontal pipes</li>
</ul>
It does not model changes in pressure resulting from significant variation of flow velocity along the flow path (with the assumption of a constant cross sectional area it must result from fluid density changes, such as in two-phase flow).
 
When connecting two components, e.g. two pipes, the momentum balance across the connection point reduces to</p> 
<pre>pipe1.port_b.p = pipe2.port_a.p</pre>
<p>This is only true if the flow velocity remains the same on each side of the connection. For any significant change in diameter (and if the resulting effects, such as change in kinetic energy, cannot be neglected) an adapter component should be used. This also allows for taking into account friction losses with respect to the actual geometry of the connection point.</p>
 
</html>",
      revisions="<html>
<ul>
<li><i>27 Nov 2008</i>
    by R&uuml;diger Franke:<br>
       Major overhaul</li>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  end DistributedPipe_Old;

  partial model PartialPipe_Old "Base class for pipe models"
    extends Modelica_Fluid.Interfaces.PartialTwoPort;

     //Initialization
    parameter Medium.AbsolutePressure p_a_start=system.p_start
      "Start value of pressure at port a" 
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_b_start=p_a_start
      "Start value of pressure at port b" 
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.MassFlowRate m_flow_start = system.m_flow_start
      "Start value for mass flow rate" 
       annotation(Evaluate=true, Dialog(tab = "Initialization"));

    //Geometry
    parameter SI.Length length "Length"   annotation(Dialog(tab="General", group="Geometry"));
    parameter SI.Diameter diameter "Diameter of circular pipe"      annotation(Dialog(group="Geometry", enable=isCircular));
    parameter SI.Length roughness(min=0)=2.5e-5
      "Average height of surface asperities (default = smooth steel pipe)" 
        annotation(Dialog(group="Geometry",enable=WallFriction.use_roughness));
    parameter Boolean isCircular=true
      "= true if cross sectional area is circular" 
      annotation (Evaluate, Dialog(tab="General", group="Geometry"));
    parameter SI.Length perimeter=Modelica.Constants.pi*diameter
      "Inner perimeter"                                                                                       annotation(Dialog(tab="General", group="Geometry", enable=not isCircular));
    parameter SI.Area crossArea=Modelica.Constants.pi*diameter*diameter/4
      "Inner cross section area"            annotation(Dialog(tab="General", group="Geometry", enable=not isCircular));
    final parameter SI.Volume V=crossArea*length "volume size";

    // Static head
    parameter SI.Length height_ab=0.0 "Height(port_b) - Height(port_a)" 
                                                                   annotation(Dialog(group="Static head"), Evaluate=true);

    // Pressure loss
    replaceable package WallFriction = 
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent
      constrainedby
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.PartialWallFriction
      "Characteristic of wall friction"  annotation(Dialog(group="Pressure loss"), choicesAllMatching=true);
     parameter SI.Diameter diameter_h=4*crossArea/perimeter
      "Hydraulic diameter"                                       annotation(Dialog(tab="General", group="Pressure loss"));
    parameter Boolean use_eta_nominal = false
      "= true, if eta_nominal is used, otherwise computed from medium" 
           annotation(Dialog(tab="Advanced", group="Pressure loss"),Evaluate=true);
    parameter SI.DynamicViscosity eta_nominal=Medium.dynamicViscosity(Medium.setState_pTX(Medium.p_default, Medium.T_default, Medium.X_default))
      "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
        annotation(Dialog(tab="Advanced", group="Pressure loss",enable=use_eta_nominal));
    parameter Boolean use_d_nominal=false
      "= true, if d_nominal is used, otherwise computed from medium"                                annotation(Dialog(tab="Advanced", group="Pressure loss"),Evaluate=true);
    parameter SI.Density d_nominal = Medium.density_pTX(Medium.p_default, Medium.T_default, Medium.X_default)
      "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
       annotation(Dialog(tab="Advanced", group="Pressure loss",enable=use_d_nominal));

    // Heat transfer
    replaceable model HeatTransfer = 
        Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_constAlpha 
      constrainedby
      Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer
      "Wall heat transfer model" 
      annotation (Dialog(group="Heat transfer"),editButton=false,choicesAllMatching=true);

    parameter Boolean show_Re = false
      "= true, if Reynolds number is included for plotting" 
       annotation (Evaluate=true, Dialog(tab="Advanced"));
    parameter Boolean from_dp=true
      " = true, use m_flow = f(dp), otherwise dp = f(m_flow)" 
      annotation (Evaluate=true, Dialog(tab="Advanced", group="Pressure loss"));
    parameter SI.AbsolutePressure dp_small = 1
      "Within regularization if |dp| < dp_small (may be wider for large discontinuities in static head)"
      annotation(Dialog(tab="Advanced", group="Pressure loss", enable=from_dp and WallFriction.use_dp_small));
    parameter SI.MassFlowRate m_flow_small = 0.01
      "Within regularizatio if |m_flow| < m_flow_small (may be wider for large discontinuities in static head)"
      annotation(Dialog(tab="Advanced", group="Pressure loss", enable=not from_dp and WallFriction.use_m_flow_small));

    annotation (defaultComponentName="pipe",Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={
          Rectangle(
            extent={{-100,44},{100,-44}},
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

  end PartialPipe_Old;

partial model PartialDistributedFlow_Old
    "Base class for a finite volume flow model"
    import Modelica_Fluid.Types;

  extends Modelica_Fluid.Vessels.BaseClasses.PartialDistributedVolume(final n = nNodes);
  extends Modelica_Fluid.ToBeRemoved.PartialPipe_Old;

//Discretization
  parameter Integer nNodes(min=1)=1 "Number of discrete flow volumes";

//Advanced model options
  parameter Boolean useInnerPortProperties=false
      "=true, port properties for pressure drop correlation are taken from neighboring control volume"
                                                                                                     annotation(Dialog(tab="Advanced", group="Pressure loss"),Evaluate=true);
//Flow quantities
  Medium.MassFlowRate[n + 1] m_flow(each min=if allowFlowReversal then -Modelica.Constants.inf else 
              0, each start=m_flow_start, each fixed=false)
      "Mass flow rates of fluid across segment boundaries";
  SI.Velocity[n+1] v "Velocity at volume boundaries (not used in balances)";
  Medium.MassFlowRate[n + 1,Medium.nXi] mXi_flow
      "Independent mass flow rates across segment boundaries";
  Medium.EnthalpyFlowRate[n + 1] H_flow
      "Enthalpy flow rates of fluid across segment boundaries";

  Medium.AbsolutePressure[n] p = medium.p "Pressure states";

  //Source terms, have to be set in inheriting class (to zero if not used)
  protected
  SI.Density[n] d=if use_d_nominal then ones(n)*d_nominal else medium.d;
  SI.Density d_a=if use_d_nominal then d_nominal else (if useInnerPortProperties then d[1] else Medium.density_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)));
  SI.Density d_b=if use_d_nominal then d_nominal else (if useInnerPortProperties then d[n] else Medium.density_phX(port_b.p, inStream(port_b.h_outflow), inStream(port_b.Xi_outflow)));
  public
  Interfaces.HeatPorts_a[nNodes] heatPorts annotation (Placement(transformation(extent={{-10,44},
            {10,64}}),                                                                               iconTransformation(extent={{-30,44},{32,60}})));
equation
  // Boundary conditions
  port_a.h_outflow = medium[1].h;
  port_b.h_outflow = medium[n].h;
  port_a.m_flow    = m_flow[1];
  port_b.m_flow    = -m_flow[n + 1];
  port_a.C_outflow = inStream(port_b.C_outflow);
  port_b.C_outflow = inStream(port_a.C_outflow);

  // Distributed flow quantities, upwind discretization
  for i in 2:n loop
    H_flow[i] = semiLinear(m_flow[i], medium[i - 1].h, medium[i].h);
    mXi_flow[i, :] = semiLinear(m_flow[i], medium[i - 1].Xi, medium[i].Xi);
    v[i] = m_flow[i]/(medium[i - 1].d + medium[i].d)*2/crossArea;
  end for;
  H_flow[1] = semiLinear(port_a.m_flow, inStream(port_a.h_outflow), medium[1].h);
  H_flow[n + 1] = -semiLinear(port_b.m_flow, inStream(port_b.h_outflow), medium[n].h);
  mXi_flow[1, :] = semiLinear(port_a.m_flow, inStream(port_a.Xi_outflow), medium[1].Xi);
  mXi_flow[n + 1, :] = -semiLinear(port_b.m_flow, inStream(port_b.Xi_outflow), medium[n].Xi);
  v[1] = m_flow[1]/(d_a + d[1])*2/crossArea;
  v[n + 1] = m_flow[n + 1]/(d[n] + d_b)*2/crossArea;

  //Mass and energy balances
  //dynamic mass balances, n "thermal" states, n pressure states (if not singleState_hydraulic)
  for i in 1:n loop
    ms_flow[i] = m_flow[i] - m_flow[i + 1];
    msXi_flow[i, :] = mXi_flow[i, :] - mXi_flow[i + 1, :];
  end for;
  //dynamic energy balances, n "thermal" states, n pressure states (if not singleState_hydraulic)
  for i in 1:n loop
    Hs_flow[i] = H_flow[i] - H_flow[i + 1];
  end for;

  annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
              -100},{100,100}}),
                       graphics),
                        Icon(coordinateSystem(preserveAspectRatio=true,
            extent={{-100,-100},{100,100}}), graphics={Rectangle(
            extent={{-100,40},{100,-40}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255})}),
      Documentation(info="<html>
<p>The model <b>PartialDistributedFlow</b> is used as a base class for pipe flows with one-dimensional spatial discretization according to the finite volume method. The flow path is divided into <tt><b>nNodes</b></tt> segments.</p>
 
<p><b>Mass and energy balances</b></p>
<p>One total mass and one energy balance is formed across each segment. If the medium contains more than one component, substance mass balances are added. Changes in potential and kinetic energy are neglected in the energy balance. The following source (or sink) terms are used in the balances and must be specified in extending models to complete this partial class:</p>
<ul>
<li>Energy balance: <tt><b>Qs_flow</b></tt>, e.g. convective or latent heat flow rate across segment boundary, and <tt><b>Ws_flow</b></tt>, e.g. mechanical power</li>
<li>Total mass balance: <tt><b>ms_flow</b></tt>, e.g. condensing mass flow of negligible volume such as water in moist air</li>
<li>Substance mass balance: <tt><b>msXi_flow</b></tt>, as above</li>
</ul>
If the <tt>dynamicsType</tt> is <b>DynamicsType.SteadyState</b> then no mass or energy is stored in the component and the mass and energy balances are reduced to a quasi steady-state formulation. It should be noted that dynamic balances are required if flow reversal should be allowed.
<p>An extending class shall define volume vector <tt><b>Vi</b></tt>, which specifies the volume of each segment, and the mass flow rates <tt><b>m_flow</b></tt> through a mementum balance.
 
<p><b>Momentum balance</b></p>
<p>The momentum balance needs to be defined by an extending class.</p> 
 
</html>", revisions="<html>
<ul>
<li><i>27 Nov 2008</i>
    by R&uuml;diger Franke:<br>
       Major overhaul</li>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
end PartialDistributedFlow_Old;

  model MassFlowRatio "simple flow multiplier"
    replaceable package Medium=Modelica.Media.Interfaces.PartialMedium annotation(choicesAllMatching);
    parameter Integer nOutlets=1
      "Number of outlet ports (mass is distributed evenly between the outlet ports";
    Modelica_Fluid.Interfaces.FluidPort_a port_a(
                                  redeclare package Medium=Medium) 
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}},
            rotation=0)));
    Modelica_Fluid.Interfaces.FluidPorts_b ports_b[nOutlets] 
                                    annotation (Placement(transformation(extent=
             {{90,-40},{110,40}}, rotation=0)));

    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics={
          Line(
            points={{-80,0},{80,0}},
            color={0,128,255},
            thickness=1),
          Line(
            points={{-80,0},{80,28}},
            color={0,128,255},
            thickness=1),
          Line(
            points={{-80,0},{80,-28}},
            color={0,128,255},
            thickness=1),
          Text(
            extent={{-150,100},{150,60}},
            lineColor={0,0,255},
            textString="%name")}),
                            Documentation(info="<html>
<p>
This model describes a simple flow partitioning, which is very helpful in cases where the flow is evenly distributed to several parallel flow paths which are identical in their dimensions and boundary conditions, as e.g. in heat exchangers. Only one of the parallel pipes needs to be simulated then. All flow variables in <b>port_b[i]</b> are equal to those at <b>port_a</b> divided by <b>nOutlets</b>. All effort variables are equal at all ports.
</p>
</html>"),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
              {100,100}}),
              graphics));

  equation
    port_a.h_outflow  =  sum( {inStream(ports_b[i].h_outflow)     for i in 1:nOutlets})/nOutlets;
    port_a.Xi_outflow = {sum( {inStream(ports_b[i].Xi_outflow[j]) for i in 1:nOutlets})/nOutlets for j in 1:Medium.nXi};
    port_a.C_outflow  = {sum( {inStream(ports_b[i].C_outflow[j])  for i in 1:nOutlets})/nOutlets for j in 1:Medium.nXi};

    for i in 1:nOutlets loop
       ports_b[i].h_outflow  = inStream(port_a.h_outflow);
       ports_b[i].Xi_outflow = inStream(port_a.Xi_outflow);
       ports_b[i].C_outflow  = inStream(port_a.C_outflow);

       // Momentum balance
       port_a.p = ports_b[i].p;

       // Mass balance
       ports_b[i].m_flow = -port_a.m_flow/nOutlets;
    end for;
  end MassFlowRatio;

  model HeatFlowRatio "simple heat flow multiplier"
    parameter Integer nOutlets=1
      "Number of outlet ports (heat is distributed evenly between the outlet ports";
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics={
          Line(
            points={{-80,0},{80,30}},
            color={127,0,0},
            thickness=1),
          Line(
            points={{-80,0},{80,0}},
            color={127,0,0},
            thickness=1),
          Line(
            points={{-80,0},{80,-32}},
            color={127,0,0},
            thickness=1),
          Text(
            extent={{-150,100},{150,60}},
            lineColor={0,0,255},
            textString="%name")}),
                            Documentation(info="<html>
<p>
Simple model for heat flow partitioning between the two ports. The heat flow rate in port_a is divided by parameter <b>nOutlets</b> to achieve the rate at ports port_b. All temperatures are equal. The model may be used e.g. for parallel pipes in heat exchangers.
</p>
</html>"));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort_a 
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}},
            rotation=0)));
    Modelica_Fluid.Interfaces.HeatPorts_b heatPorts_b[nOutlets] 
                                   annotation (Placement(transformation(extent=
              {{90,-40},{110,40}}, rotation=0)));
  equation
    for i in 1:nOutlets loop
       heatPorts_b[i].Q_flow = -heatPort_a.Q_flow/nOutlets;
       heatPorts_b[i].T      =  heatPort_a.T;
    end for;
  end HeatFlowRatio;

model StaticHead
    "Models the static head between two ports at different heights"
  extends Modelica_Fluid.PressureLosses.BaseClasses.PartialTwoPortTransport;
  parameter SI.Length height_ab "Height(port_b) - Height(port_a)";
  parameter Medium.MassFlowRate m_flow_small(min=0) = 1e-4
      "For bi-directional flow, density is regularized in the region |m_flow| < m_flow_small (m_flow_small > 0 required)"
    annotation(Dialog(tab="Advanced"));

  Medium.Density d "Density of the passing fluid";
  outer Modelica_Fluid.System system "System properties";
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics={
          Rectangle(
            extent={{-100,60},{100,-60}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={192,192,192}),
          Rectangle(
            extent={{-100,48},{100,-48}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255}),
          Text(
            extent={{-150,80},{150,120}},
            lineColor={0,0,255},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255},
            textString="%name")}),         Documentation(info="<html>
<p>
This model describes the static head due to the relative height between the two connectors. No mass, energy and momentum storage, and no pressure drop due to friction are considered.
</p>
</html>",
        revisions="<html>
<ul>
<li><i>31 Oct 2007</i>
    by <a href=\"mailto:jonas@modelon.se\">Jonas Eborn</a>:<br>
       Changed to flow-direction dependent density</li>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Added to Modelica_Fluid</li>
</ul>
</html>"));
equation
//  d = if dp > 0 then medium_a.d else medium_b.d;
  /* Density is currently not handled correctly. The correct formula:
        d = if port_a.m_flow>0 then d_a else d_b
     leads to a non-linear system of equations that mays have an infinite number of
     solutions and also no solution may exist for small flow rates.
     Most likely, this can only be correctly handled with a dynamic momentum balance.
  */
  if allowFlowReversal then
     d = (port_a_d_inflow + port_b_d_inflow)/2
        "temporary solution that must be improved (BUT NOT WITH if port_a.m_flow>0 then d_a else d_b)";
  else
     d = port_a_d_inflow;
  end if;
  dp = height_ab*system.g*d;
end StaticHead;

  model JunctionVolume
    "Splitting/joining component with static balances for a dynamic control volume"
    import Modelica_Fluid.Types;
    import Modelica_Fluid.Types.PortFlowDirection;
    outer Modelica_Fluid.System system "System properties";

    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
      "Fluid medium model" 
        annotation (choicesAllMatching=true);
    parameter SI.Volume V "Volume";

    // Assumptions
    parameter Modelica_Fluid.Types.Dynamics dynamicsType=system.dynamicsType
      "Dynamics option" 
      annotation(Evaluate=true, Dialog(tab = "Assumptions"));

    SI.InternalEnergy U "Internal energy";
    SI.Mass m "Total mass";
    SI.Mass[Medium.nXi] mXi "Independent masses";

    Modelica_Fluid.Interfaces.FluidPort_a port_1(
      redeclare package Medium=Medium,
      m_flow(min=if (portFlowDirection_1==PortFlowDirection.Entering) then 0.0 else -Modelica.Constants.inf,
      max=if (portFlowDirection_1==PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) 
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}},
            rotation=0)));
    Modelica_Fluid.Interfaces.FluidPort_b port_2(
      redeclare package Medium=Medium,
      m_flow(min=if (portFlowDirection_2==PortFlowDirection.Entering) then 0.0 else -Modelica.Constants.inf,
      max=if (portFlowDirection_2==PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) 
      annotation (Placement(transformation(extent={{90,-10},{110,10}}, rotation=
             0)));
    Modelica_Fluid.Interfaces.FluidPort_a port_3(
      redeclare package Medium=Medium,
      m_flow(min=if (portFlowDirection_3==PortFlowDirection.Entering) then 0.0 else -Modelica.Constants.inf,
      max=if (portFlowDirection_3==PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) 
      annotation (Placement(transformation(extent={{-10,90},{10,110}}, rotation=
             0)));

    Medium.ExtraProperty C[Medium.nC] "Trace substance mixture content";
    Medium.BaseProperties medium(preferredMediumStates=true);

    parameter Types.Init initType=system.initType "Initialization option" 
      annotation(Evaluate=true,Dialog(tab="Initialization"));
    parameter Medium.AbsolutePressure p_start "Start value of pressure" 
      annotation(Dialog(tab="Initialization"));
    parameter Boolean use_T_start=true "=true, use T_start, otherwise h_start" 
      annotation(Dialog(tab="Initialization"),Evaluate=true);
    parameter Medium.Temperature T_start=
      if use_T_start then system.T_start else Medium.temperature_phX(p_start,h_start,X_start)
      "Start value of temperature" 
      annotation(Dialog(tab="Initialization",enable=use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.specificEnthalpy_pTX(p_start,T_start,X_start) else Medium.h_default
      "Start value of specific enthalpy" 
      annotation(Dialog(tab="Initialization",enable=not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX]=Medium.X_default
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization",enable=Medium.nXi>0));

  protected
    parameter PortFlowDirection portFlowDirection_1=PortFlowDirection.Bidirectional
      "Flow direction for port_1" 
     annotation(Dialog(tab="Advanced"));
    parameter PortFlowDirection portFlowDirection_2=PortFlowDirection.Bidirectional
      "Flow direction for port_2" 
     annotation(Dialog(tab="Advanced"));
    parameter PortFlowDirection portFlowDirection_3=PortFlowDirection.Bidirectional
      "Flow direction for port_3" 
     annotation(Dialog(tab="Advanced"));

    annotation (Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={
          Rectangle(
            extent={{-100,41},{100,-47}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={192,192,192}),
          Rectangle(
            extent={{-100,37},{100,-43}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255}),
          Rectangle(
            extent={{-34,100},{34,37}},
            lineColor={0,0,0},
            fillPattern=FillPattern.VerticalCylinder,
            fillColor={192,192,192}),
          Rectangle(
            extent={{-30,100},{30,35}},
            lineColor={0,0,0},
            fillPattern=FillPattern.VerticalCylinder,
            fillColor={0,127,255}),
          Ellipse(
            extent={{-9,10},{11,-10}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-150,-60},{150,-100}},
            lineColor={0,0,255},
            textString="%name")}),
      Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics));
  initial equation
    // Initial conditions
    if initType == Types.Init.NoInit then
      // no initial equations
    elseif initType == Types.Init.InitialValues then
      medium.p = p_start;
      medium.h = h_start;
    elseif initType == Types.Init.SteadyState then
      der(medium.p) = 0;
      der(medium.h) = 0;
    elseif initType == Types.Init.SteadyStateHydraulic then
      der(medium.p) = 0;
      medium.h = h_start;
    else
      assert(false, "Unsupported initialization option");
    end if;

  equation
    // Only one connection allowed to a port to avoid unwanted ideal mixing
  /*
  assert(cardinality(port_1) <= 1,"
port_1 of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller.
");
  assert(cardinality(port_2) <= 1,"
port_2 of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller.
");
  assert(cardinality(port_3) <= 1,"
port_3 of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller.
");
*/

    // Boundary conditions
    port_1.h_outflow = medium.h;
    port_2.h_outflow = medium.h;
    port_3.h_outflow = medium.h;

    port_1.Xi_outflow = medium.Xi;
    port_2.Xi_outflow = medium.Xi;
    port_3.Xi_outflow = medium.Xi;

    // Internal quantities
    m   = medium.d*V;
    mXi = m*medium.Xi;
    U   = m*medium.u;

    // Mass balances
    if dynamicsType < Types.Dynamics.SteadyStateMass then
      der(m)   = port_1.m_flow + port_2.m_flow + port_3.m_flow "Mass balance";
      der(mXi) = port_1.m_flow*actualStream(port_1.Xi_outflow)
                + port_2.m_flow*actualStream(port_2.Xi_outflow)
                + port_3.m_flow*actualStream(port_3.Xi_outflow)
        "Component mass balances";
    else
         0   = port_1.m_flow + port_2.m_flow + port_3.m_flow "Mass balance";
      zeros(Medium.nXi) = port_1.m_flow*actualStream(port_1.Xi_outflow)
                + port_2.m_flow*actualStream(port_2.Xi_outflow)
                + port_3.m_flow*actualStream(port_3.Xi_outflow)
        "Component mass balances";
    end if;

  /* 
  zeros(Medium.nC) = port_1.m_flow*actualStream(port_1.C_outflow)
                      + port_2.m_flow*actualStream(port_2.C_outflow)
                      + port_3.m_flow*actualStream(port_3.C_outflow) 
    "Trace substance mass balances";
*/

    // Momentum balance (suitable for compressible media)
    port_1.p = medium.p;
    port_2.p = medium.p;
    port_3.p = medium.p;

    // Energy balance
    if dynamicsType < Types.Dynamics.SteadyState then
      der(U) = port_1.m_flow*actualStream(port_1.h_outflow)
              + port_2.m_flow*actualStream(port_2.h_outflow)
              + port_3.m_flow*actualStream(port_3.h_outflow);
    else
        0  = port_1.m_flow*actualStream(port_1.h_outflow)
              + port_2.m_flow*actualStream(port_2.h_outflow)
              + port_3.m_flow*actualStream(port_3.h_outflow);
    end if;
  end JunctionVolume;

  model GenericJunction
    "Branching component with balances for a dynamic control volume"
    import Modelica.Constants;
    import Modelica_Fluid.Types;
    import Modelica_Fluid.Types.ModelStructure;
    outer Modelica_Fluid.System system "System properties";
    parameter Integer nPorts_a(min=1)=1 "Number of ports on side a";
    parameter Integer nPorts_b(min=1)=1 "Number of ports on side b";
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
      "Fluid medium model" 
        annotation (choicesAllMatching=true);
    parameter SI.Volume V "Volume";
    parameter SI.Pressure dp_nominal "nominal (linear) pressure drop" annotation(Dialog(enable=not modelStructure==ModelStructure.av_vb));
    parameter SI.MassFlowRate m_flow_nominal "nominal mass flow rate"  annotation(Dialog(enable=not modelStructure==ModelStructure.av_vb));

    SI.InternalEnergy U "Internal energy";
    SI.Mass m "Total mass";
    SI.Mass[Medium.nXi] mXi "Independent masses";

    Interfaces.FluidPorts_a[nPorts_a] ports_a(
      redeclare each package Medium=Medium,
      m_flow(each min=if allowFlowReversal then -Constants.inf else 0))
      "Fluid connectors a (positive design flow direction is from ports_a to ports_b)"
      annotation (Placement(
          transformation(extent={{-110,40},{-90,-40}}, rotation=0)));
    Interfaces.FluidPorts_b[nPorts_b] ports_b(
      redeclare each package Medium=Medium,
      m_flow(each max=if allowFlowReversal then +Constants.inf else 0))
      "Fluid connectors b (positive design flow direction is from ports_a to ports_b)"
      annotation (Placement(
          transformation(extent={{90,40},{110,-40}}, rotation=0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort "Thermal port"
      annotation (Placement(transformation(extent={{-10,90},{10,110}}, rotation=0)));
    Medium.ExtraProperty C[Medium.nC] "Trace substance mixture content";
    Medium.BaseProperties medium(T(start=T_start),p(start=p_start),h(start=h_start),X(start=X_start), preferredMediumStates=true);

    // Assumptions
    parameter Boolean allowFlowReversal = system.allowFlowReversal
      "allow flow reversal, false restricts to design direction (port_a -> port_b)"
      annotation(Dialog(tab="Assumptions"), Evaluate=true);
    parameter Modelica_Fluid.Types.Dynamics dynamicsType=system.dynamicsType
      "Dynamics option" 
      annotation(Evaluate=true, Dialog(tab = "Assumptions"));

    // Initialization
    parameter Types.Init initType=Types.Init.NoInit "Initialization option" 
      annotation(Evaluate=true,Dialog(tab="Initialization"));
    parameter Medium.AbsolutePressure p_start "Start value of pressure" 
      annotation(Dialog(tab="Initialization"));
    parameter Boolean use_T_start=true "=true, use T_start, otherwise h_start" 
      annotation(Dialog(tab="Initialization"),Evaluate=true);
    parameter Medium.Temperature T_start=
      if use_T_start then system.T_start else Medium.temperature_phX(p_start,h_start,X_start)
      "Start value of temperature" 
      annotation(Dialog(tab="Initialization",enable=use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.specificEnthalpy_pTX(p_start,T_start,X_start) else Medium.h_default
      "Start value of specific enthalpy" 
      annotation(Dialog(tab="Initialization",enable=not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX]=Medium.X_default
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization",enable=Medium.nXi>0));

    parameter ModelStructure modelStructure=ModelStructure.av_vb annotation(Evaluate=true);

    Medium.EnthalpyFlowRate ports_a_H_flow[nPorts_a];
    Medium.EnthalpyFlowRate ports_b_H_flow[nPorts_b];
    Medium.MassFlowRate ports_a_mXi_flow[nPorts_a,Medium.nXi];
    Medium.MassFlowRate ports_b_mXi_flow[nPorts_b,Medium.nXi];
    Medium.ExtraPropertyFlowRate ports_a_mC_flow[nPorts_a,Medium.nC];
    Medium.ExtraPropertyFlowRate ports_b_mC_flow[nPorts_b,Medium.nC];

    annotation (Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={
          Ellipse(
            extent={{-19,0},{1,-20}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillPattern=FillPattern.Sphere,
            fillColor={0,128,255}),
          Ellipse(
            extent={{-9,10},{11,-10}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-150,150},{150,110}},
            lineColor={0,0,255},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255},
            textString="%name")}),
      Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics));

  initial equation
    // Initial conditions
    if initType == Types.Init.NoInit then
      // no initial equations
    elseif initType == Types.Init.InitialValues then
      medium.p = p_start;
      medium.h = h_start;
    elseif initType == Types.Init.SteadyState then
      der(medium.p) = 0;
      der(medium.h) = 0;
    elseif initType == Types.Init.SteadyStateHydraulic then
      der(medium.p) = 0;
      medium.h = h_start;
    else
      assert(false, "Unsupported initialization option");
    end if;

  equation
    // Only one connection allowed to a port to avoid unwanted ideal mixing
  /*
for i in 1:nPorts_a loop
  assert(cardinality(ports_a[i]) <= 1,"
ports_a[" + String(i) + "] of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller.
");
end for;
 
for i in 1:nPorts_b loop
  assert(cardinality(ports_b[i]) <= 1,"
ports_a[" + String(i) + "] of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller.
");
end for;
*/

    heatPort.T = medium.T;

    if dynamicsType < Types.Dynamics.SteadyStateMass then
      sum(ports_a.m_flow)+sum(ports_b.m_flow) = der(m) "Mass balance";

      for i in 1:Medium.nXi loop
        sum(ports_a_mXi_flow[:,i])+sum(ports_b_mXi_flow[:,i]) = der(mXi[i])
          "Substance mass balance";
      end for;
    else
      sum(ports_a.m_flow)+sum(ports_b.m_flow) = 0 "Mass balance";

      for i in 1:Medium.nXi loop
        sum(ports_a_mXi_flow[:,i])+sum(ports_b_mXi_flow[:,i]) = 0
          "Substance mass balance";
      end for;
    end if;

    if dynamicsType < Types.Dynamics.SteadyState then
      sum(ports_a_H_flow) + sum(ports_b_H_flow) + heatPort.Q_flow = der(U)
        "Energy balance";
    else
      sum(ports_a_H_flow) + sum(ports_b_H_flow) + heatPort.Q_flow = 0
        "Energy balance";
    end if;

    for i in 1:Medium.nC loop
      sum(ports_a_mC_flow[:,i])+sum(ports_b_mC_flow[:,i]) = 0
        "Trace substance mass balance";
    end for;

    for i in 1:nPorts_a loop
      ports_a[i].h_outflow  = medium.h;
      ports_a[i].Xi_outflow = medium.Xi;
      ports_a[i].C_outflow = C;

      ports_a_H_flow[i] = ports_a[i].m_flow * actualStream(ports_a[i].h_outflow)
        "Enthalpy flow";
      ports_a_mXi_flow[i,:] = ports_a[i].m_flow * actualStream(ports_a[i].Xi_outflow)
        "Component mass flow";
      ports_a_mC_flow[i,:] = ports_a[i].m_flow * actualStream(ports_a[i].C_outflow)
        "Trace substance mass flow";
    end for;

    for i in 1:nPorts_b loop
      ports_b[i].h_outflow  = medium.h;
      ports_b[i].Xi_outflow = medium.Xi;
      ports_b[i].C_outflow = C;

      ports_b_H_flow[i] = ports_b[i].m_flow * actualStream(ports_b[i].h_outflow)
        "Enthalpy flow";
      ports_b_mXi_flow[i,:] = ports_b[i].m_flow * actualStream(ports_b[i].Xi_outflow)
        "Component mass flow";
      ports_b_mC_flow[i,:] = ports_b[i].m_flow * actualStream(ports_b[i].C_outflow)
        "Trace substance mass flow";
    end for;

    if modelStructure==ModelStructure.av_vb or modelStructure == ModelStructure.av_b then
      ports_a.p=fill(medium.p, nPorts_a);
    else
      ports_a.p-fill(medium.p,nPorts_a) = ports_a.m_flow*dp_nominal/m_flow_nominal;
    end if;

    if modelStructure==ModelStructure.av_vb or modelStructure==ModelStructure.a_vb then
      ports_b.p=fill(medium.p,nPorts_b);
    else
      ports_b.p-fill(medium.p,nPorts_b)=ports_b.m_flow*dp_nominal/m_flow_nominal;
    end if;

    U=m*medium.u;
    mXi=m*medium.Xi;
    m=medium.d*V;

  end GenericJunction;

  model Test1 "Prescribed inputs, initial values"

    Modelica_Fluid.Examples.HeatingSystem plant(
        tank(initType=Modelica_Fluid.Types.Init.InitialValues),
        pipe(initType=Modelica_Fluid.Types.Init.InitialValues),
        radiator(initType=Modelica_Fluid.Types.Init.InitialValues)) 
                annotation (Placement(transformation(extent={{0,0},{20,20}},
            rotation=0)));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}),
                        graphics),
      experiment(StopTime=6000, Tolerance=1e-006),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Initial equations with initial values for the states close to the steady state are selected for all components.
<p>
The simulation initializes and runs for 6000 seconds without problems.
</html>"));
  end Test1;

  model Test2 "Prescribed inputs, all derivatives equal to zero"
    extends Modelica_Fluid.ToBeRemoved.Test1(
                  plant(
        tank(initType=Modelica_Fluid.Types.Init.SteadyState),
        pipe(initType=Modelica_Fluid.Types.Init.SteadyState),
        radiator(initType=Modelica_Fluid.Types.Init.SteadyState)));
    annotation (
      Documentation(info="<html>
Initial equations for steady-state are selected for all components.
<p>
The simulation initializes and runs. However, due to the mathematical structure of the model, the initial level and temperature of the tank are not defined: any initial level and temperature in fact would satisfy the initial equations. Therefore, it is not possible for the user to actually set those values precisely. Accordingly, it turns out that when a simulation is run, the initial values of those two states in the simulation depend in an unpredictable fashion from the start values; in other words, if the start value for the level and temperature are changed, the initial values also change: this should not happen if the steady-state initialization problem had a well-defined, unique solution, at least if the start values are chosen close to the solution. 
</html>
"),   experiment(StopTime=6000, Tolerance=1e-006),
      experimentSetupOutput(equdistant=false));
  end Test2;

  model Test3 "Prescribed inputs, all derivatives equal to zero"
    extends Modelica_Fluid.ToBeRemoved.Test1(
                  plant(
        tank(initType=Modelica_Fluid.Types.Init.SteadyState),
        pipe(initType=Modelica_Fluid.Types.Init.SteadyState),
        radiator(initType=Modelica_Fluid.Types.Init.SteadyState)));
    annotation (
      Documentation(info="<html>
Initial equations for steady-state are selected for all components. Moreover, additional initial equations are given to specify initial values for the tank temperature and tank pressure are specified at the system level (this is necessary because the tank components currently does not support that initialization option).
<p>The resulting initialization problem has more equations than unknowns states, but still has a unique solution. Unfortunately, Dymola 6.0d cannot find it. Also note that the equation Dymola suggests to remove is in fact necessary to determine the plant.valve.Av parameter.
</html>
"),   experiment(StopTime=6000, Tolerance=1e-006),
      experimentSetupOutput(equdistant=false));
  initial equation
   // plant.tank.level = 1;
   // plant.tank.medium.T = Modelica.SIunits.Conversions.from_degC(20);
  end Test3;

  model Test4
    "Prescribed inputs, all derivatives equal to zero for the pipes, initial states given for the tank"
    extends Modelica_Fluid.ToBeRemoved.Test1(
                  plant(
        tank(initType=Modelica_Fluid.Types.Init.InitialValues),
        pipe(initType=Modelica_Fluid.Types.Init.SteadyState),
        radiator(initType=Modelica_Fluid.Types.Init.SteadyState)));
    annotation (
      Documentation(info="<html>
Initial equations for steady-state are selected for the pipe components, initial values are specified for the tank. Due to the mathematical structure of the tank model, steady state conditions can be obtained with arbitrary initial pressure and temperature, so this initial condition equations actually correspond to a steady-state condition. However, a user which does not know the exact mathematical structure of the model could not be aware of that.
<p>
</html>
"),   experiment(StopTime=6000, Tolerance=1e-006),
      experimentSetupOutput(equdistant=false));
  end Test4;

  model Test5
    "Prescribed inputs, all derivatives equal to zero, zero pressure loss in the radiator"
    extends Modelica_Fluid.ToBeRemoved.Test1(
                  plant(
        tank(initType=Modelica_Fluid.Types.Init.SteadyState),
        pipe(initType=Modelica_Fluid.Types.Init.SteadyState),
        radiator(initType=Modelica_Fluid.Types.Init.SteadyState,
                 redeclare model PressureDrop = 
              Modelica_Fluid.Pipes.BaseClasses.PressureDrop.NominalPressureDrop(dp_nominal=0,smoothFlowReversal=true))));

    annotation (
      Documentation(info="<html>
Initial equations for steady-state are selected for all components. The model of the radiator pipe has zero pressure losses.
<p>
The radiator pipe has no pressure losses in the momentum balances, so the pressure corresponding to its PortVolume component is equal to the pressure of the tank, leading to an index 2 DAE. One pressure state is removed by the dummy derivative algorithm, then Dymola complains there is one initial equation too many. In fact, the initialization problem is overdetermined, but still has infinitely many solutions. Also note that the equation Dymola suggests to remove is in fact necessary to determine the plant.valve.Av parameter.
</html>
"),   experiment(StopTime=6000, Tolerance=1e-006),
      experimentSetupOutput(equdistant=false));
  end Test5;

  model Test6
    "Prescribed inputs, all derivatives equal to zero, zero pressure loss in the radiator"
    extends Modelica_Fluid.ToBeRemoved.Test1(
                  plant(
        tank(initType=Modelica_Fluid.Types.Init.SteadyState),
        pipe(initType=Modelica_Fluid.Types.Init.SteadyState),
        radiator(initType=Modelica_Fluid.Types.Init.SteadyState,
                 redeclare model PressureDrop = 
              Modelica_Fluid.Pipes.BaseClasses.PressureDrop.NominalPressureDrop(dp_nominal=0,smoothFlowReversal=true))));
    annotation (
      Documentation(info="<html>
Initial equations for steady-state are selected for all components, plus additional initial equations to set the initial level and temperature of the tank. The model of the radiator pipe has zero pressure losses.
<p>
The radiator pipe has no pressure losses in the momentum balances, so the pressure corresponding to its PortVolume component is equal to the pressure of the tank, leading to an index 2 DAE. One pressure state is removed by the dummy derivative algorithm, then Dymola complains there are three initial equations too many. In fact, the initialization problem is overdetermined, but has one unique solution. Also note that the equation Dymola suggests to remove is in fact necessary to determine the plant.valve.Av parameter.
</html>
"),   experiment(StopTime=6000, Tolerance=1e-006),
      experimentSetupOutput(equdistant=false));
  initial equation
    plant.tank.level = 1;
    plant.tank.medium.T = Modelica.SIunits.Conversions.from_degC(20);
  end Test6;

  model Test7
    "Prescribed inputs, carefully selected initial conditions, zero pressure loss in the radiator"
    extends Modelica_Fluid.ToBeRemoved.Test1(
                  plant(
        tank(initType=Modelica_Fluid.Types.Init.InitialValues),
        pipe(initType=Modelica_Fluid.Types.Init.SteadyState),
        radiator(redeclare model PressureDrop = 
              Modelica_Fluid.Pipes.BaseClasses.PressureDrop.NominalPressureDrop(dp_nominal=0,smoothFlowReversal=true),
            initType=Modelica_Fluid.Types.Init.NoInit)));
    annotation (
      Documentation(info="<html>
Initial equations for steady-state are selected for all components, plus additional initial equations to set the initial level and temperature of the tank. The model of the radiator pipe has zero pressure losses.
<p>
The initial conditions here have been carefully crafted: 
<ul>
<li>Steady-state conditions for the tank have been replaced by initial states, because otherwise the problem has multiple solutions
<li>Steady-state conditions for the radiator have been replaced by a steady-state condition only on the thermal state (enthalpy), not on the pressure that has been removed by the dummy derivative algorithm. 
</ul>
The simulation initializes and runs correctly.
<p>
However, it is apparent that such a special treatment is out of the question for more complex models, whose mathematical structure might be much harder to understand. It is also extremely inconvenient that the user has to change initial conditions depending on index reduction: a user shouldn't need to know the inner mathematical details of a model (in this case, the fact that there is actually a state at the port of two directly connected components, due to a special selection of the WallFriction parameter) in order to be able to run it.
<p> 
On the other hand, it is quite straightforward to set up overdetermined conditions such as those of Test6, corresponding to steadyState for every model, plus additional equations to actually set the tank level and temperature. The latter can be easily added when one discovers that the steadyState option is not enough to get the desired initial values of those two variables.
</html>
"),   experiment(StopTime=6000, Tolerance=1e-006),
      experimentSetupOutput(equdistant=false));
  initial equation
    der(plant.radiator.volume.medium.h) = 0;
  end Test7;
end ToBeRemoved;
