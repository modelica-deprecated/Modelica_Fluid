within Modelica_Fluid;
package Pipes "Lumped, distributed and thermal pipe components"
    extends Modelica_Fluid.Icons.VariantLibrary;

  model StaticPipe "Basic pipe flow model without storage of mass or energy"
    extends Modelica_Fluid.Pipes.BaseClasses.PartialPipe(
                                           redeclare model HeatTransfer = 
          BaseClasses.HeatTransfer.PipeHT_ideal,                                                               final
        area_h =                                                                                                    0);
    PressureLosses.WallFrictionAndGravity wallFriction(
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
  end StaticPipe;

  model LumpedPipe
    // Assumptions
    parameter Modelica_Fluid.Types.Dynamics dynamicsType=system.dynamicsType
      "Dynamics option" 
      annotation(Evaluate=true, Dialog(tab = "Assumptions"));

    // Initialization
    parameter Types.Init initType=system.initType "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));

    // Extend here to get right ordering in parameter box
    extends Modelica_Fluid.Pipes.BaseClasses.PartialPipe;

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

    HeatTransfer heatTransfer(
      redeclare final package Medium = Medium,
      final diameter_h=diameter_h,
      final area_h=area_h,
      final crossArea=crossArea,
      final length=length,
      state=volume.medium.state,
      m_flow = 0.5*(port_a.m_flow - port_b.m_flow),
      final useFluidHeatPort=true) 
      annotation (Placement(transformation(extent={{-20,0},{20,40}},   rotation=0)));

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
    Modelica_Fluid.Volumes.Volume volume(
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
    connect(heatPort, heatTransfer.wallHeatPort) annotation (Line(
        points={{0,54},{0,34}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(heatTransfer.fluidHeatPort, volume.heatPort) annotation (Line(
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
  end LumpedPipe;

 model DistributedPipeLumpedPressure
    "Distributed pipe model with lumped pressure state"
   import Modelica_Fluid.Types.ModelStructure;
   extends Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow(
     Qs_flow=heatTransfer.Q_flow,
     Ws_flow=zeros(n),
     ms_flow=zeros(n),
     msXi_flow=zeros(n, Medium.nXi),
     Vi=ones(n)*V/n,
     final port_a_exposesState = (modelStructure == ModelStructure.av_b) or (modelStructure == ModelStructure.avb),
     final port_b_exposesState = (modelStructure == ModelStructure.a_vb) or (modelStructure == ModelStructure.avb));

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
   HeatTransfer[n] heatTransfer(
     redeclare each final package Medium = Medium,
     each final diameter_h=diameter_h,
     each final area_h=area_h/n,
     each final crossArea=crossArea,
     each final length=length/n,
     state=medium.state,
     m_flow = 0.5*(m_flow[1:n]+m_flow[2:n+1])) "Convective heat transfer" 
             annotation (Placement(transformation(extent={{-20,-5},{20,35}},  rotation=0)));

   SI.Length[2] dlength "discretized length for pressure drop";
   SI.Length[2] dheight_ab "discretized height_ab for static head";

  protected
   SI.DynamicViscosity eta_a=if not WallFriction.use_eta then 1.e-10 else (if 
       use_eta_nominal then eta_nominal else (if use_approxPortProperties then Medium.dynamicViscosity(medium[1].state) else Medium.dynamicViscosity(Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)))));
   SI.DynamicViscosity eta_b=if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then eta_nominal else (if use_approxPortProperties then Medium.dynamicViscosity(medium[n].state) else Medium.dynamicViscosity(Medium.setState_phX(port_b.p, inStream(port_b.h_outflow), inStream(port_b.Xi_outflow)))));

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
   else // avb
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
     Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}},
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

 end DistributedPipeLumpedPressure;

  model DistributedPipe "Distributed pipe model"
  import Modelica_Fluid.Types.ModelStructure;
  extends Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow(
    Qs_flow=heatTransfer.Q_flow,
    Ws_flow=zeros(n),
    ms_flow=zeros(n),
    msXi_flow=zeros(n, Medium.nXi),
    Vi=ones(n)*V/n,
    final port_a_exposesState = (modelStructure == ModelStructure.av_b) or (modelStructure == ModelStructure.avb),
    final port_b_exposesState = (modelStructure == ModelStructure.a_vb) or (modelStructure == ModelStructure.avb));

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
  HeatTransfer[n] heatTransfer(
    redeclare each final package Medium = Medium,
    each final diameter_h=diameter_h,
    each final area_h=area_h/n,
    each final crossArea=crossArea,
    each final length=length/n,
    state=medium.state,
    m_flow = 0.5*(m_flow[1:n]+m_flow[2:n+1])) "Convective heat transfer" 
              annotation (Placement(transformation(extent={{-20,-5},{20,35}},  rotation=0)));

  SI.Length[n+1] dlength "discretized length for pressure drop";
  SI.Length[n+1] dheight_ab "discretized height_ab for static head";
  protected
  SI.DynamicViscosity eta_a=if not WallFriction.use_eta then 1.e-10 else (if 
      use_eta_nominal then eta_nominal else (if use_approxPortProperties then eta[1] else Medium.dynamicViscosity(Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)))));
  SI.DynamicViscosity eta_b=if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then eta_nominal else (if use_approxPortProperties then eta[n] else Medium.dynamicViscosity(Medium.setState_phX(port_b.p, inStream(port_b.h_outflow), inStream(port_b.Xi_outflow)))));
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

    //Pressure drop and gravity
    //Simplified Momentum Balance, dp contains contributions from gravitational and friction effects
    dp = cat(1, {port_a.p}, p) - cat(1, p, {port_b.p});

    if modelStructure == ModelStructure.a_v_b then
      dlength = cat(1, {0.5*length/n}, length/n*ones(n-1), {0.5*length/n});
    elseif modelStructure == ModelStructure.av_b then
      dlength = cat(1, {0}, length/n*ones(n));
    elseif modelStructure == ModelStructure.a_vb then
      dlength = cat(1, length/n*ones(n), {0});
    else // avb
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
<li><tt>avb</tt>: n-1 momentum balances between first and nth volume, potential pressure states at both ports. It's use should be avoided, since not the entire pipe length is taken into account.
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
  end DistributedPipe;

  package BaseClasses
    extends Modelica_Fluid.Icons.BaseClassLibrary;

    partial model PartialPipe "Base class for one dimensional flow models"
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
        "Convective heat transfer" 
        annotation (Dialog(group="Heat transfer"),editButton=true,choicesAllMatching=true);
      parameter SI.Area area_h=perimeter*length "Heat transfer area" 
                                                               annotation(Dialog(tab="General", group="Heat transfer"));

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

    end PartialPipe;

  partial model PartialDistributedFlow
      "Base class for a finite volume flow model"
      import Modelica_Fluid.Types;

    parameter Modelica_Fluid.Types.Dynamics dynamicsType=system.dynamicsType
        "Dynamics option" 
      annotation(Evaluate=true, Dialog(tab = "Assumptions"));

  //Discretization
    parameter Integer nNodes(min=1)=1 "Number of discrete flow volumes";
    final parameter Integer n = nNodes;

    final parameter Boolean static = dynamicsType == Types.Dynamics.SteadyState
        "= true, static balances, no mass or energy is stored" 
                                  annotation(Dialog(tab="Assumptions"),Evaluate=true);

  //Initialization
    parameter Types.Init initType=system.initType "Initialization option" 
        annotation(Evaluate=true, Dialog(tab = "Initialization"));

    // Extend here to get right ordering in parameter box
    extends Modelica_Fluid.Pipes.BaseClasses.PartialPipe;

    final parameter Medium.AbsolutePressure[n] p_start=if n > 1 then linspace(
          p_a_start - (p_a_start - p_b_start)/(2*n),
          p_b_start + (p_a_start - p_b_start)/(2*n),
          n) else {(p_a_start + p_b_start)/2} "Start value of pressure";

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

  //Advanced model options
    parameter Boolean use_approxPortProperties=false
        "=true, port properties for pressure drop correlation are taken from neighboring control volume"
                                                                                                       annotation(Dialog(tab="Advanced", group="Pressure loss"),Evaluate=true);
    input SI.Volume[n] Vi "Discretized volume, determine in inheriting class ";

  //Total quantities
    SI.Energy[n] U "Internal energy of fluid";
    SI.Mass[n] m "Fluid mass";
    SI.Mass[n,Medium.nXi] mXi "Substance mass";

  //Flow quantities
    Medium.MassFlowRate[n + 1] m_flow(each min=if allowFlowReversal then -Modelica.Constants.inf else 
                0, each start=m_flow_start, each fixed=false)
        "Mass flow rates of fluid across segment boundaries";
    SI.Velocity[n+1] v "Velocity at volume boundaries (not used in balances)";
    Medium.MassFlowRate[n + 1,Medium.nXi] mXi_flow
        "Independent mass flow rates across segment boundaries";
    Medium.EnthalpyFlowRate[n + 1] H_flow
        "Enthalpy flow rates of fluid across segment boundaries";

    Medium.BaseProperties[n] medium(
      each preferredMediumStates=if static then false else true,
      p(start=p_start),
      each h(start=h_start),
      each T(start=T_start),
      each Xi(start=X_start[1:Medium.nXi]));
    Medium.AbsolutePressure[n] p = medium.p "Pressure states";

    //Source terms, have to be set in inheriting class (to zero if not used)
    protected
    input Medium.MassFlowRate[n] ms_flow "Mass flow rate, source or sink";
    input Medium.MassFlowRate[n,Medium.nXi] msXi_flow
        "Independent mass flow rates, source or sink";
    input SI.HeatFlowRate[n] Qs_flow "Heat flow rate, source or sink";
    input SI.Power[n] Ws_flow "Mechanical power, p*der(V) etc.";
    SI.Density[n] d=if use_d_nominal then ones(n)*d_nominal else medium.d;
    SI.Density d_a=if use_d_nominal then d_nominal else (if use_approxPortProperties then d[1] else Medium.density_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)));
    SI.Density d_b=if use_d_nominal then d_nominal else (if use_approxPortProperties then d[n] else Medium.density_phX(port_b.p, inStream(port_b.h_outflow), inStream(port_b.Xi_outflow)));
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

    // Total quantities
    for i in 1:n loop
      m[i] =Vi[i]*medium[i].d;
      mXi[i, :] = m[i]*medium[i].Xi;
      U[i] = m[i]*medium[i].u;
    end for;

    //Mass and energy balances
    if dynamicsType < Types.Dynamics.SteadyStateMass then
    //dynamic mass balances, n "thermal" states, n pressure states (if not singleState_hydraulic)
      for i in 1:n loop
        der(m[i]) = m_flow[i] - m_flow[i + 1] + ms_flow[i];
        der(mXi[i, :]) = mXi_flow[i, :] - mXi_flow[i + 1, :] + msXi_flow[i, :];
      end for;
    else
    //steady state mass balances, no numerical states, no flow reversal possible
      for i in 1:n loop
        0 = m_flow[i] - m_flow[i + 1] + ms_flow[i];
        zeros(Medium.nXi) = mXi_flow[i, :] - mXi_flow[i + 1, :] + msXi_flow[i, :];
      end for;
    end if;
    if dynamicsType < Types.Dynamics.SteadyState then
    //dynamic energy balances, n "thermal" states, n pressure states (if not singleState_hydraulic)
      for i in 1:n loop
        der(U[i]) = H_flow[i] - H_flow[i + 1] + Qs_flow[i];
      end for;
    else
    //steady state energy balances, no numerical states, no flow reversal possible
      for i in 1:n loop
        0 = H_flow[i] - H_flow[i + 1] + Qs_flow[i];
      end for;
    end if;
    for i in 1:n loop
      assert((allowFlowReversal and not static) or (m_flow[i] >= 0), "Flow reversal not allowed in distributed pipe");
    end for;

  initial equation
    // Initial conditions
    if not static then
      if initType == Types.Init.NoInit then
      // no initial equations
      elseif initType == Types.Init.SteadyState then
      //steady state initialization
        if use_T_start then
          der(medium.T) = zeros(n);
        else
          der(medium.h) = zeros(n);
        end if;
        if not (Medium.singleState) then
          der(medium.p) = zeros(n);
        end if;
        for i in 1:n loop
          der(medium[i].Xi) = zeros(Medium.nXi);
        end for;
      elseif initType == Types.Init.InitialValues then
      //Initialization with initial values
        if use_T_start then
          medium.T = ones(n)*T_start;
        else
          medium.h = ones(n)*h_start;
        end if;
        if not Medium.singleState then
           medium.p=p_start;
        end if;
      elseif initType == Types.Init.SteadyStateHydraulic then
      //Steady state initialization for hydraulic states (p)
        if use_T_start then
          medium.T = ones(n)*T_start;
        else
          medium.h = ones(n)*h_start;
        end if;
        if not Medium.singleState then
          der(medium.p) = zeros(n);
        end if;
      else
        assert(false, "Unsupported initialization option");
      end if;
    end if;

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
 
</html>",   revisions="<html>
<ul>
<li><i>27 Nov 2008</i>
    by R&uuml;diger Franke:<br>
       Major overhaul</li>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  end PartialDistributedFlow;

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

  package HeatTransfer
    partial model PartialPipeHeatTransfer
        "base class for any pipe heat transfer correlation"

      // heat port to wall
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a wallHeatPort
          "Thermal port to wall" 
        annotation (Placement(transformation(extent={{-10,60},{10,80}},
                rotation=0)));

      // Parameters
      replaceable package Medium=Modelica.Media.Interfaces.PartialMedium annotation(Dialog(tab="No input", enable=false));
      parameter SI.Area area_h "Total heat transfer area" annotation(Dialog(tab="No input", enable=false));
      parameter SI.Length diameter_h "Hydraulic diameter" annotation(Dialog(tab="No input", enable=false));
      parameter SI.Area crossArea "Cross flow area" annotation(Dialog(tab="No input", enable=false));
      parameter SI.Length length "Pipe length" annotation(Dialog(tab="No input", enable=false));

      // Inputs provided to heat transfer model
      input Medium.ThermodynamicState state;
      input SI.MassFlowRate m_flow;

      // Output defined by heat transfer model
      output SI.HeatFlowRate Q_flow "Heat flow rates";

      // Internal variables
      SI.Temperature T;
      parameter Boolean useFluidHeatPort = false
          "= true to use fluidHeatPort instead of output Q_flow" 
        annotation(Dialog(tab="No input", enable=false), Evaluate=true, HideResult=true);
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b fluidHeatPort if useFluidHeatPort
          "Thermal port to fluid" 
        annotation (Placement(transformation(extent={{-10,-70},{10,-50}},
                rotation=0)));
      Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow
          prescribedHeatFlow "Needed to connect to conditional connector" 
        annotation (Placement(transformation(extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={0,-10})));

    equation
      T = Medium.temperature(state);
      if useFluidHeatPort then
        prescribedHeatFlow.Q_flow = Q_flow;
      else
        prescribedHeatFlow.port.T = T;
      end if;
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
      connect(prescribedHeatFlow.port, fluidHeatPort) annotation (Line(
          points={{-1.83697e-015,-20},{0,-20},{0,-60}},
          color={191,0,0},
          smooth=Smooth.None));
    end PartialPipeHeatTransfer;

    partial model PartialPipeHT_Nu
        "Base class for pipe heat transfer correlation in terms of Nusselt numberheat transfer in a circular pipe for laminar and turbulent one-phase flow"
      extends PartialPipeHeatTransfer;
      parameter SI.CoefficientOfHeatTransfer alpha0=100;
      SI.CoefficientOfHeatTransfer alpha(each start=alpha0)
          "CoefficientOfHeatTransfer";
      Real Re "Reynolds number";
      Real Pr "Prandtl number";
      Real Nu "Nusselt number";
      SI.DynamicViscosity eta "Dynamic viscosity";
      SI.ThermalConductivity lambda "Thermal conductivity";
    equation
      eta=Medium.dynamicViscosity(state);
      lambda=Medium.thermalConductivity(state);
      Pr = Medium.prandtlNumber(state);
      Re = CharacteristicNumbers.ReynoldsNumber(m_flow, diameter_h, crossArea, eta);
      Nu = CharacteristicNumbers.NusseltNumber(alpha, diameter_h, lambda);
      wallHeatPort.Q_flow=Q_flow;
      wallHeatPort.Q_flow=alpha*area_h*(wallHeatPort.T - T);
        annotation (Documentation(info="<html>
Base class for heat transfer models that are expressed in terms of the Nusselt number and which can be used in distributed pipe models.
</html>"));
    end PartialPipeHT_Nu;

    model PipeHT_ideal "Ideal heat transfer without thermal resistance"
      extends PartialPipeHeatTransfer;
    equation
      T = wallHeatPort.T;
      Q_flow = wallHeatPort.Q_flow;
      annotation(Documentation(info="<html>
Ideal heat transfer without thermal resistance.
</html>"));
    end PipeHT_ideal;

    model PipeHT_constAlpha
      extends PartialPipeHeatTransfer;
      parameter SI.CoefficientOfHeatTransfer alpha0=200;
      annotation(Documentation(info="<html>
Simple heat transfer correlation with constant heat transfer coefficient, used as default component in <a distributed pipe models.
</html>"));
    equation
      wallHeatPort.Q_flow = alpha0*area_h*(wallHeatPort.T - T);
      wallHeatPort.Q_flow = Q_flow;
    end PipeHT_constAlpha;
    annotation (Documentation(info="<html>
Heat transfer correlations for pipe models
</html>"));

    model PipeHT_LamTurb_local
        "laminar and turbulent forced convection in pipes, local coeff."
      extends PartialPipeHT_Nu;
      protected
      Real Nu_turb "Nusselt number for turbulent flow";
      Real Nu_lam "Nusselt number for laminar flow";
      Real Nu_1;
      Real Nu_2;
      Real Xi;
    equation
      Nu_1=3.66;
      Nu_turb=smooth(0,(Xi/8)*abs(Re)*Pr/(1+12.7*(Xi/8)^0.5*(Pr^(2/3)-1))*(1+1/3*(diameter_h/length)^(2/3)));
      Xi=(1.8*Modelica.Math.log10(max(1e-10,Re))-1.5)^(-2);
      Nu_lam=(Nu_1^3+0.7^3+(Nu_2-0.7)^3)^(1/3);
      Nu_2=smooth(0,1.077*(abs(Re)*Pr*diameter_h/length)^(1/3));
      Nu=spliceFunction(Nu_turb, Nu_lam, Re-6150, 3850);
      annotation (Documentation(info="<html>
Heat transfer model for laminar and turbulent flow in pipes. Range of validity:
<ul>
<li>fully developed pipe flow</li>
<li>forced convection</li>
<li>one phase Newtonian fluid</li>
<li>(spatial) constant wall temperature in the laminar region</li>
<li>0 &le; Re &le; 1e6, 0.6 &le; Pr &le; 100, d/L &le; 1</li>
<li>The correlation holds for non-circular pipes only in the turbulent region. Use diameter_h=4*area/perimeter as characteristic length.</li>
</ul>
The correlation takes into account the spatial position along the pipe flow, which changes discontinuously at flow reversal. However, the heat transfer coefficient itself is continuous around zero flow rate, but not its derivative.
<h4><font color=\"#008000\">References</font></h4>
 
<dl><dt>Verein Deutscher Ingenieure (1997):</dt>
    <dd><b>VDI W&auml;rmeatlas</b>.
         Springer Verlag, Ed. 8, 1997.</dd>
</dl>
</html>"));
    end PipeHT_LamTurb_local;
  end HeatTransfer;

  end BaseClasses;
  annotation (Documentation(info="<html>
 
</html>"));

end Pipes;
