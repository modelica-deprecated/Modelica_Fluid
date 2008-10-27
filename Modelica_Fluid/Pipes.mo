within Modelica_Fluid;
package Pipes "Lumped, distributed and thermal pipe components"
    extends Modelica_Fluid.Icons.VariantLibrary;

  model StaticPipe "Basic pipe flow model without storage of mass or energy"
    extends PressureLosses.WallFrictionAndGravity;
  end StaticPipe;

  model LumpedPipe
    "Short pipe with one volume, wall friction and gravity effect"
    outer Modelica_Fluid.System system "System properties";
    replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
      annotation (choicesAllMatching = true);
    parameter Modelica_Fluid.Types.Init initType=Modelica_Fluid.Types.Init.NoInit
      "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_a_start
      "Start value of pressure at port_a" 
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_b_start = p_a_start
      "Start value of pressure at port_b" 
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.MassFlowRate m_flow_start = 0
      "Guess value of port_a.m_flow" 
      annotation(Dialog(tab = "Initialization"));
    parameter Boolean use_T_start = true
      "= true, use T_start, otherwise h_start" 
      annotation(Dialog(tab = "Initialization"), Evaluate=true);
    parameter Medium.Temperature T_start=
      if use_T_start then Medium.T_default else Medium.temperature_phX(p_a_start,h_start,X_start)
      "Start value of temperature" 
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.specificEnthalpy_pTX(p_a_start, T_start, X_start) else Medium.h_default
      "Start value of specific enthalpy" 
      annotation(Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    replaceable package WallFriction = 
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent
      constrainedby
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.PartialWallFriction
      "Characteristic of wall friction"  annotation(choicesAllMatching=true);

    parameter Modelica.SIunits.Length length "Length of pipe";
    parameter Modelica.SIunits.Diameter diameter
      "Inner (hydraulic) diameter of pipe";
    parameter Modelica.SIunits.Length height_ab=0.0
      "Height(port_b) - Height(port_a)"                                annotation(Evaluate=true);
    parameter Modelica.SIunits.Length roughness(min=0)=2.5e-5
      "Absolute roughness of pipe (default = smooth steel pipe)" 
        annotation(Dialog(enable=WallFriction.use_roughness));
    parameter Boolean use_nominal = false
      "= true, if eta_nominal and d_nominal are used, otherwise computed from medium"
        annotation(Evaluate=true);
    parameter Modelica.SIunits.DynamicViscosity eta_nominal=0.01
      "Nominal dynamic viscosity (for wall friction computation)" annotation(Dialog(enable=use_nominal));
    parameter Modelica.SIunits.Density d_nominal=0.01
      "Nominal density (for wall friction computation)" annotation(Dialog(enable=use_nominal));
    parameter Modelica_Fluid.Types.FlowDirection flowDirection=
        system.flowDirection
      "Unidirectional (port_a -> port_b) or bidirectional flow component" 
       annotation(Dialog(tab="Advanced"));
    parameter Modelica.SIunits.AbsolutePressure dp_small=1
      "Turbulent flow if |dp| >= dp_small (only used if WallFriction=QuadraticTurbulent)"
      annotation(Dialog(tab="Advanced", enable=WallFriction.use_dp_small));

      Modelica_Fluid.Interfaces.FluidPort_a port_a(
                                    redeclare package Medium = Medium)
      "Fluid connector a (positive design flow direction is from port_a to port_b)"
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}},
            rotation=0)));
      Modelica_Fluid.Interfaces.FluidPort_b port_b(
                                    redeclare package Medium = Medium)
      "Fluid connector b (positive design flow direction is from port_a to port_b)"
        annotation (Placement(transformation(extent={{110,-10},{90,10}},
            rotation=0)));

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
            extent={{-150,-94},{150,-134}},
            lineColor={0,0,255},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255},
            textString="%name"),
          Ellipse(
            extent={{-11,10},{9,-10}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Line(
            points={{40,-70},{-50,-70}},
            color={0,128,255},
            smooth=Smooth.None),
          Polygon(
            points={{30,-55},{70,-70},{30,-85},{30,-55}},
            lineColor={0,128,255},
            smooth=Smooth.None,
            fillColor={0,128,255},
            fillPattern=FillPattern.Solid)}),Documentation(info="<html>
<p>
Simple pipe model consisting of one volume, 
wall friction (with different friction correlations)
and gravity effect. This model is mostly used to demonstrate how
to build up more detailed models from the basic components.
Note, if the \"thermalPort\" is not connected, then the pipe
is totally insulated (= no thermal flow from the fluid to the
pipe wall/environment).
</p>
</html>"),
      Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics),
      uses(Modelica_Fluid(version="1.0 Beta 2"), Modelica(version="2.2.2")));
    Modelica_Fluid.PressureLosses.WallFrictionAndGravity wallFriction1(
      redeclare package Medium = Medium,
      flowDirection=flowDirection,
      redeclare package WallFriction = WallFriction,
      length=length/2,
      diameter=diameter,
      height_ab=height_ab/2,
      roughness=roughness,
      use_nominal=use_nominal,
      eta_nominal=eta_nominal,
      d_nominal=d_nominal,
      from_dp=true,
      dp_small=dp_small,
      show_Re=false,
      dp_start = (p_a_start - p_b_start)/2,
      m_flow_start = m_flow_start,
      compute_T=false) 
      annotation (Placement(transformation(extent={{-60,-10},{-40,10}},
            rotation=0)));
    Volumes.MixingVolume volume(
      redeclare package Medium = Medium,
      V=Modelica.Constants.pi*(diameter/2)^2*length,
      flowDirection=flowDirection,
      initType=initType,
      p_start=(p_a_start+p_b_start)/2,
      use_T_start=use_T_start,
      T_start=T_start,
      h_start=h_start,
      X_start=X_start) 
      annotation (Placement(transformation(extent={{-10,-10},{10,10}}, rotation=
             0)));
    Modelica_Fluid.PressureLosses.WallFrictionAndGravity wallFriction2(
      redeclare package Medium = Medium,
      flowDirection=flowDirection,
      redeclare package WallFriction = WallFriction,
      length=length/2,
      diameter=diameter,
      height_ab=height_ab/2,
      roughness=roughness,
      use_nominal=use_nominal,
      eta_nominal=eta_nominal,
      d_nominal=d_nominal,
      from_dp=true,
      dp_small=dp_small,
      show_Re=false,
      dp_start = (p_a_start - p_b_start)/2,
      m_flow_start = m_flow_start,
      compute_T=false)              annotation (Placement(transformation(extent=
             {{40,-10},{60,10}}, rotation=0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a thermalPort 
      annotation (Placement(transformation(extent={{-10,44},{10,64}}, rotation=
              0)));
  equation
    connect(wallFriction1.port_a, port_a) 
      annotation (Line(points={{-60,0},{-100,0}}, color={0,127,255}));
    connect(wallFriction2.port_b, port_b) 
      annotation (Line(points={{60,0},{100,0}}, color={0,127,255}));
    connect(volume.thermalPort, thermalPort) 
      annotation (Line(points={{0,9.8},{0,54}}, color={191,0,0}));
    connect(volume.port_a, wallFriction1.port_b)       annotation (Line(
        points={{-10.2,0},{-40,0}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(volume.port_b, wallFriction2.port_a)       annotation (Line(
        points={{10,0},{40,0}},
        color={0,127,255},
        smooth=Smooth.None));
  end LumpedPipe;

 model DistributedPipeLumpedPressure "Distributed pipe model"

   extends
      Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlowLumpedPressure(
     Qs_flow=heatTransfer.Q_flow,
     Ws_flow=zeros(n),
     ms_flow=zeros(n),
     msXi_flow=zeros(n, Medium.nXi),
     Vi=ones(n)*V/n);

    //Volume size
   final parameter SI.Volume V=area*length;
 //Heat transfer area
   parameter SI.Area area_h=perimeter*length "Heat transfer area" 
                                                              annotation(Dialog(tab="General", group="Heat transfer"));
    //Pressure Drop
   replaceable package WallFriction = 
       Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent
     constrainedby
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.PartialWallFriction
      "Characteristic of wall friction" annotation(Dialog(tab="General", group="Pressure loss"),choicesAllMatching=true);
   parameter SI.Diameter d_h=4*area/perimeter "Hydraulic diameter" 
                                                                  annotation(Dialog(tab="General", group="Pressure loss"));
   parameter SI.Length height_ab=0.0 "Height(port_b) - Height(port_a)" 
                                                                     annotation(Dialog(tab="General", group="Geometry"),Evaluate=true);
   parameter SI.Length roughness(min=0) = 2.5e-5
      "Absolute roughness of pipe (default = smooth steel pipe)" 
      annotation(Dialog(tab="General", group="Pressure loss",enable=WallFriction.use_roughness));
   parameter Boolean use_eta_nominal=false
      "= true, if eta_nominal is used, otherwise computed from medium"                            annotation(Dialog(tab="General", group="Pressure loss"),Evaluate=true);

   parameter SI.DynamicViscosity eta_nominal=0.01
      "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
                                                                            annotation(Dialog(tab="General", group="Pressure loss",enable=use_nominal));
   parameter Boolean show_Re=false
      "= true, if Reynolds number is included for plotting" 
     annotation (Evaluate=true, Dialog(tab="Advanced", group="Pressure loss"));
   parameter Boolean from_dp=true
      " = true, use m_flow = f(dp), otherwise dp = f(m_flow)" 
    annotation (Evaluate=true, Dialog(tab="Advanced", group="Pressure loss"));
   parameter SI.AbsolutePressure dp_small=1
      "Turbulent flow if |dp| >= dp_small (only used if WallFriction=QuadraticTurbulent)"
    annotation(Dialog(tab="Advanced",group="Pressure loss", enable=from_dp and WallFriction.use_dp_small));

   parameter SI.MassFlowRate m_flow_small=0.01
      "Turbulent flow if |m_flow| >= m_flow_small (only used if WallFriction=QuadraticTurbulent)"
    annotation(Dialog(tab="Advanced",group="Pressure loss", enable=not from_dp and WallFriction.use_m_flow_small));

   SI.ReynoldsNumber[n+1] Re=Modelica_Fluid.Utilities.ReynoldsNumber_m_flow(
       m_flow,
       (cat(1, {eta_a}, eta) + cat(1, eta, {eta_b}))*0.5,
       diameter) if                                                                                              show_Re
      "Reynolds number of pipe flow";
   inner Medium.ThermodynamicState[n] state;
   replaceable Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_constAlpha
      heatTransfer(
     redeclare final package Medium = Medium,
     final n=n,
     final d_h=d_h,
     final A_h=area_h,
     final A_cross=area,
     final L=length,
     T=medium.T) constrainedby
      Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer(
     redeclare package Medium = Medium,
     n=n,
     d_h=d_h,
     A_h=area_h,
     A_cross=area,
     L=length,
     T=medium.T) "Convective heat transfer" 
             annotation (Dialog(tab="General", group="Heat transfer"),editButton=true,choicesAllMatching,
      Placement(transformation(extent={{-20,-20},{20,20}}, rotation=0)));
   Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort
      "Thermal port" 
 annotation (Placement(transformation(extent={{-10,44},{10,64}}, rotation=0)));
   outer Modelica_Fluid.System system "System properties";

  protected
   SI.DynamicViscosity eta_a=if not WallFriction.use_eta then 1.e-10 else (if 
       use_eta_nominal then eta_nominal else (if use_approxPortProperties then Medium.dynamicViscosity(medium[1].state) else Medium.dynamicViscosity(Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)))));
   SI.DynamicViscosity eta_b=if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then eta_nominal else (if use_approxPortProperties then Medium.dynamicViscosity(medium[n].state) else Medium.dynamicViscosity(Medium.setState_phX(port_b.p, inStream(port_b.h_outflow), inStream(port_b.Xi_outflow)))));

   annotation (
     Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}},
          grid={1,1}), graphics={
          Rectangle(
            extent={{-100,44},{100,40}},
            lineColor={0,0,0},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-100,-40},{100,-44}},
            lineColor={0,0,0},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
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
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-150,-97},{150,-137}},
            lineColor={0,0,255},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255},
            textString="%name"),
          Line(
            points={{30,-70},{-60,-70}},
            color={0,128,255},
            smooth=Smooth.None),
          Polygon(
            points={{20,-55},{60,-70},{20,-85},{20,-55}},
            lineColor={0,128,255},
            smooth=Smooth.None,
            fillColor={0,128,255},
            fillPattern=FillPattern.Solid)}),
     Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}},
          grid={1,1}),
             graphics),
     Documentation(info="<html>
Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow_pLumped\">PartialDistributedFlow_pLumped</a>. Source terms in mass and energy balances are set to zero. The total volume is a paramter. The number of momentum balances is reduced to two, one on each side of the hydraulic state, which corresponds to a constant pressure along the entire pipe with pressure drop and gravitational forces lumped at the ports.<The additional component <tt>heatTransfer</tt> specifies the source term <tt>Qs_flow</tt> in the energy balance. The default component uses a constant coefficient of heat transfer to model convective heat transfer between segment boundary (<tt>thermalPort</tt>) and the bulk flow. The <tt>heatTransfer</tt> model is replaceable and can be exchanged with any model extended from <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer\">PartialPipeHeatTransfer</a>. .
</html>",  revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));

 equation
   state=medium.state;

 if from_dp and not WallFriction.dp_is_zero then
     m_flow[1] = WallFriction.massFlowRate_dp(
       dp[1] - ((integer(n/2) + 1)*2 - 1)/(2*n)*height_ab*system.g*(d[1]+d[n])
         /2,
       d_a,
       d_b,
       eta_a,
       eta_b,
       ((integer(n/2) + 1)*2 - 1)/(2*n)*length,
       d_h,
       roughness,
       dp_small);
     m_flow[n + 1] = WallFriction.massFlowRate_dp(
       dp[2] - (2*n - (integer(n/2) + 1)*2 + 1)/(2*n)*height_ab*system.g*(d[1]+d[n])/2,
       d_a,
       d_b,
       eta_a,
       eta_b,
       (2*n - (integer(n/2) + 1)*2 + 1)/(2*n)*length,
       d_h,
       roughness,
       dp_small);
 else
   dp[1] = (WallFriction.pressureLoss_m_flow(
       m_flow[1],
       d_a,
       d_b,
       eta_a,
       eta_b,
       ((integer(n/2) + 1)*2 - 1)/(2*n)*length,
       d_h,
       roughness,
       m_flow_small) + ((integer(n/2) + 1)*2 - 1)/(2*n)*height_ab*system.g*(
       d[1]+d[n])/2);
     dp[2] = (WallFriction.pressureLoss_m_flow(
       m_flow[n+1],
       d_a,
       d_b,
       eta_a,
       eta_b,
       (2*n - (integer(n/2) + 1)*2 + 1)/(2*n)*length,
       d_h,
       roughness,
       m_flow_small) + (2*n - (integer(n/2) + 1)*2 + 1)/(2*n)*height_ab*
       system.g*(d[1]+d[n])/2);
 end if;

   connect(thermalPort, heatTransfer.thermalPort) 
     annotation (Line(points={{0,54},{0,14}}, color={191,0,0}));
 end DistributedPipeLumpedPressure;

  model DistributedPipe "Distributed pipe model"
  extends Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow(
    Qs_flow=heatTransfer.Q_flow,
    Ws_flow=zeros(n),
    ms_flow=zeros(n),
    msXi_flow=zeros(n, Medium.nXi),
    Vi=ones(n)*V/n,
    final modelStructure=Modelica_Fluid.Types.ModelStructure.a_v_b,
    redeclare final Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare
          package Medium = 
        Medium, m_flow(min=if allowFlowReversal and not static then -Modelica.Constants.inf else  0)),
    redeclare final Modelica_Fluid.Interfaces.FluidPort_b port_b(redeclare
          package Medium = 
        Medium, m_flow(max=if allowFlowReversal and not static then +Modelica.Constants.inf else 0)));

  parameter SI.Area area_h=perimeter*length "Heat transfer area" 
                                                             annotation(Dialog(tab="General", group="Heat transfer"));
  final parameter SI.Volume V=area*length;
  parameter SI.Length height_ab=0.0 "Height(port_b) - Height(port_a)" 
                                                                   annotation(Dialog(tab="General", group="Geometry"),Evaluate=true);

    //Pressure Drop
  replaceable package WallFriction = 
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent
                                                                                constrainedby
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.PartialWallFriction
      "Characteristic of wall friction" 
                                       annotation(Dialog(tab="General", group="Pressure loss"),choicesAllMatching=true);
  parameter Boolean use_d_nominal=false
      "= true, if d_nominal is used, otherwise computed from medium"                              annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);
  parameter SI.Diameter d_h=4*area/perimeter "Hydraulic diameter" 
                                                                 annotation(Dialog(tab="General", group="Pressure loss"));
  parameter SI.Length roughness(min=0) = 2.5e-5
      "Absolute roughness of pipe (default = smooth steel pipe)" 
    annotation(Dialog(tab="General", group="Pressure loss",enable=WallFriction.use_roughness));
  parameter Boolean use_eta_nominal=false
      "= true, if eta_nominal is used, otherwise computed from medium"                          annotation(Dialog(tab="General", group="Pressure loss"),Evaluate=true);
  parameter SI.DynamicViscosity eta_nominal=Medium.dynamicViscosity(Medium.setState_pTX(Medium.p_default, Medium.T_default, Medium.X_default))
      "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
                                                                          annotation(Dialog(tab="General", group="Pressure loss",enable=use_nominal));
  parameter Boolean show_Re=false
      "= true, if Reynolds number is included for plotting" 
   annotation (Evaluate=true, Dialog(tab="Advanced", group="Pressure loss"));
  parameter Boolean from_dp=true
      " = true, use m_flow = f(dp), otherwise dp = f(m_flow)" 
  annotation (Evaluate=true, Dialog(tab="Advanced", group="Pressure loss"));
  parameter SI.AbsolutePressure dp_small=1
      "Turbulent flow if |dp| >= dp_small (only used if WallFriction=QuadraticTurbulent)"
  annotation(Dialog(tab="Advanced",group="Pressure loss", enable=from_dp and WallFriction.use_dp_small));

  parameter SI.MassFlowRate m_flow_small=0.01
      "Turbulent flow if |m_flow| >= m_flow_small (only used if WallFriction=QuadraticTurbulent)"
  annotation(Dialog(tab="Advanced",group="Pressure loss", enable=not from_dp and WallFriction.use_m_flow_small));

  SI.ReynoldsNumber[n+1] Re=Modelica_Fluid.Utilities.ReynoldsNumber_m_flow(
    m_flow,
    (cat(1, {eta_a}, eta) + cat(1, eta, {eta_b}))*0.5,
    diameter) if                                                                                               show_Re
      "Reynolds number of pipe flow";
  inner Medium.ThermodynamicState[n] state;
  replaceable Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_constAlpha
      heatTransfer(
    redeclare final package Medium = Medium,
    final n=n,
    final d_h=d_h,
    final A_h=area_h,
    final A_cross=area,
    final L=length,
    T=medium.T) constrainedby
      Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer(
    redeclare package Medium = Medium,
    n=n,
    d_h=d_h,
    A_h=area_h,
    A_cross=area,
    L=length,
    T=medium.T) "Convective heat transfer" 
            annotation (Dialog(tab="General", group="Heat transfer"),editButton=true,choicesAllMatching,
      Placement(transformation(extent={{-20,-20},{20,20}}, rotation=0)));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort
      "Thermal port" annotation (Placement(transformation(extent={{-10,44},{10,
              64}}, rotation=0)));
  outer Modelica_Fluid.System system "System properties";

  protected
  SI.Pressure[n+1] dp_stat;
  SI.DynamicViscosity eta_a=if not WallFriction.use_eta then 1.e-10 else (if 
      use_eta_nominal then eta_nominal else (if use_approxPortProperties then eta[1] else Medium.dynamicViscosity(Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)))));
  SI.DynamicViscosity eta_b=if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then eta_nominal else (if use_approxPortProperties then eta[n] else Medium.dynamicViscosity(Medium.setState_phX(port_b.p, inStream(port_b.h_outflow), inStream(port_b.Xi_outflow)))));
  SI.DynamicViscosity[n] eta=if not WallFriction.use_eta then fill(1.e-10, n) else (if use_eta_nominal then fill(eta_nominal, n) else 
      Medium.dynamicViscosity(medium.state));

    annotation (
  Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={
          Rectangle(
            extent={{-100,44},{100,40}},
            lineColor={0,0,0},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-100,-40},{100,-44}},
            lineColor={0,0,0},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
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
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-150,-100},{150,-140}},
            lineColor={0,0,255},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255},
            textString="%name"),
          Line(
            points={{30,-75},{-60,-75}},
            color={0,128,255},
            smooth=Smooth.None),
          Polygon(
            points={{20,-60},{60,-75},{20,-90},{20,-60}},
            lineColor={0,128,255},
            smooth=Smooth.None,
            fillColor={0,128,255},
            fillPattern=FillPattern.Solid)}),
  Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}},
          grid={1,1}),
          graphics),
  Documentation(info="<html>
<p>Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow\">PartialDistributedFlow</a>. Source terms in the mass balances are set to zero. The total volume is a parameter. The additional component <tt>heatTransfer</tt> specifies the source term <tt>Qs_flow</tt> in the energy balance. The default component uses a constant coefficient of heat transfer to model convective heat transfer between segment boundary (<tt>thermalPort</tt>) and the bulk flow. The <tt>heatTransfer</tt> model is replaceable and can be exchanged with any model extended from <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer\">PartialPipeHeatTransfer</a>.</p>
<p>Pressure drop correlations (algebraic and possibly non-linear flow model) correlate the pressure in the first control volume with the pressure in port_a and the pressures of port_b and the nth control volume, respectively.</p>
 
</html>",
      revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));

  equation
    //Pressure difference in momentum balance
    dp = dp_stat;
  state=medium.state;

  //Pressure drop and gravity
    if from_dp and not WallFriction.dp_is_zero then
    m_flow[1] = WallFriction.massFlowRate_dp(
      dp_stat[1] - height_ab/2/n*system.g*d[1],
      d_a,
      d[1],
      eta_a,
      eta[1],
      length/n/2,
      d_h,
      roughness,
      dp_small);
    for i in 2:n loop
      m_flow[i] = WallFriction.massFlowRate_dp(
        dp_stat[i] - height_ab/n*system.g*(d[i-1] + d[i])/2,
        d[i - 1],
        d[i],
        eta[i - 1],
        eta[i],
        length/n,
        d_h,
        roughness,
        dp_small);
    end for;
    m_flow[n + 1] = WallFriction.massFlowRate_dp(
      dp_stat[n+1] - height_ab/n/2*system.g*d[n],
      d[n],
      d_b,
      eta[n],
      eta_b,
      length/n/2,
      d_h,
      roughness,
      dp_small);

    else
     dp_stat[1] = WallFriction.pressureLoss_m_flow(
        m_flow[1],
        d_a,
        d[1],
        eta_a,
        eta[1],
        length/n/2,
        d_h,
        roughness,
        m_flow_small) + height_ab/n*system.g*d[1]/2;
    for i in 2:n loop
      dp_stat[i] = WallFriction.pressureLoss_m_flow(
        m_flow[i],
        d[i - 1],
        d[i],
        eta[i - 1],
        eta[i],
        length/n,
        d_h,
        roughness,
        m_flow_small) + height_ab/n*system.g*(d[i - 1] + d[i])/2;
    end for;
    dp_stat[n] = WallFriction.pressureLoss_m_flow(
      m_flow[n+1],
      d[n],
      d_b,
      eta[n],
      eta_b,
      length/n/2,
      d_h,
      roughness,
      m_flow_small) + height_ab/n/2*system.g*d[n];

    end if;
    connect(thermalPort, heatTransfer.thermalPort) 
  annotation (Line(points={{0,54},{0,14}}, color={191,0,0}));
  end DistributedPipe;

model DistributedPipeSb "Distributed pipe model"

  extends Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow(
    Qs_flow=heatTransfer.Q_flow,
    Ws_flow=zeros(n),
    ms_flow=zeros(n),
    msXi_flow=zeros(n, Medium.nXi),
    Vi=ones(n)*V/n,
    final modelStructure=Modelica_Fluid.Types.ModelStructure.a_vb,
    redeclare final Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare
          package Medium = 
        Medium, m_flow(min=if allowFlowReversal and not static then -Modelica.Constants.inf else  0)),
    redeclare final Modelica_Fluid.Interfaces.FluidStatePort_b port_b(redeclare
          package Medium = 
        Medium, m_flow(max=if allowFlowReversal and not static then +Modelica.Constants.inf else 0)));
  parameter SI.Area area_h=perimeter*length "Heat transfer area" 
                                                             annotation(Dialog(tab="General", group="Heat transfer"));
  final parameter SI.Volume V=area*length;
  parameter SI.Length height_ab=0.0 "Height(port_b) - Height(port_a)" 
                                                                   annotation(Dialog(tab="General", group="Geometry"),Evaluate=true);

//Pressure Drop
 replaceable package WallFriction = 
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent
                                                                                constrainedby
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.PartialWallFriction
      "Characteristic of wall friction" 
                                       annotation(Dialog(tab="General", group="Pressure loss"),choicesAllMatching=true);
  parameter Boolean use_d_nominal=false
      "= true, if d_nominal is used, otherwise computed from medium"                              annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);
  parameter SI.Diameter d_h=4*area/perimeter "Hydraulic diameter" 
                                                                 annotation(Dialog(tab="General", group="Pressure loss"));
  parameter SI.Length roughness(min=0) = 2.5e-5
      "Absolute roughness of pipe (default = smooth steel pipe)" 
    annotation(Dialog(tab="General", group="Pressure loss",enable=WallFriction.use_roughness));
  parameter Boolean use_eta_nominal=false
      "= true, if eta_nominal is used, otherwise computed from medium"                          annotation(Dialog(tab="General", group="Pressure loss"),Evaluate=true);
  parameter SI.DynamicViscosity eta_nominal=Medium.dynamicViscosity(Medium.setState_pTX(Medium.p_default, Medium.T_default, Medium.X_default))
      "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
                                                                          annotation(Dialog(tab="General", group="Pressure loss",enable=use_nominal));
  parameter Boolean show_Re=false
      "= true, if Reynolds number is included for plotting" 
   annotation (Evaluate=true, Dialog(tab="Advanced", group="Pressure loss"));
  parameter Boolean from_dp=true
      " = true, use m_flow = f(dp), otherwise dp = f(m_flow)" 
  annotation (Evaluate=true, Dialog(tab="Advanced", group="Pressure loss"));
  parameter SI.AbsolutePressure dp_small=1
      "Turbulent flow if |dp| >= dp_small (only used if WallFriction=QuadraticTurbulent)"
  annotation(Dialog(tab="Advanced",group="Pressure loss", enable=from_dp and WallFriction.use_dp_small));

  parameter SI.MassFlowRate m_flow_small=0.01
      "Turbulent flow if |m_flow| >= m_flow_small (only used if WallFriction=QuadraticTurbulent)"
  annotation(Dialog(tab="Advanced",group="Pressure loss", enable=not from_dp and WallFriction.use_m_flow_small));

  SI.ReynoldsNumber[n+1] Re=Modelica_Fluid.Utilities.ReynoldsNumber_m_flow(
    m_flow,
    (cat(1, {eta_a}, eta) + cat(1, eta, {eta_b}))*0.5,
    diameter) if                                                                                               show_Re
      "Reynolds number of pipe flow";
  inner Medium.ThermodynamicState[n] state;
  replaceable Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_constAlpha
      heatTransfer(
    redeclare final package Medium = Medium,
    final n=n,
    final d_h=d_h,
    final A_h=area_h,
    final A_cross=area,
    final L=length,
    T=medium.T) constrainedby
      Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer(
    redeclare package Medium = Medium,
    n=n,
    d_h=d_h,
    A_h=area_h,
    A_cross=area,
    L=length,
    T=medium.T) "Convective heat transfer" 
            annotation (Dialog(tab="General", group="Heat transfer"),editButton=true,choicesAllMatching,
      Placement(transformation(extent={{-20,-20},{20,20}}, rotation=0)));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort
      "Thermal port" 
annotation (Placement(transformation(extent={{-10,44},{10,64}}, rotation=0)));
  outer Modelica_Fluid.System system "System properties";

  protected
  SI.Pressure[n] dp_stat;
   SI.DynamicViscosity eta_a=if not WallFriction.use_eta then 1.e-10 else (if 
      use_eta_nominal then eta_nominal else (if use_approxPortProperties then eta[1] else Medium.dynamicViscosity(Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)))));
  SI.DynamicViscosity eta_b=if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then eta_nominal else (if use_approxPortProperties then eta[n] else Medium.dynamicViscosity(Medium.setState_phX(port_b.p, inStream(port_b.h_outflow), inStream(port_b.Xi_outflow)))));
 SI.DynamicViscosity[n] eta=if not WallFriction.use_eta then fill(1.e-10, n) else (if use_eta_nominal then fill(eta_nominal, n) else 
      Medium.dynamicViscosity(medium.state));

annotation (
  Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={
          Rectangle(
            extent={{-100,44},{100,40}},
            lineColor={0,0,0},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-100,-40},{100,-44}},
            lineColor={0,0,0},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
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
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-150,-96},{150,-136}},
            lineColor={0,0,255},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255},
            textString="%name"),
          Line(
            points={{30,-75},{-60,-75}},
            color={0,128,255},
            smooth=Smooth.None),
          Polygon(
            points={{20,-60},{60,-75},{20,-90},{20,-60}},
            lineColor={0,128,255},
            smooth=Smooth.None,
            fillColor={0,128,255},
            fillPattern=FillPattern.Solid)}),
  Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}},
          grid={1,1}),
          graphics),
  Documentation(info="<html>
Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow\">PartialDistributedFlow</a>. Source terms in the mass balances are set to zero. The total volume is a parameter. The additional component <tt>heatTransfer</tt> specifies the source term <tt>Qs_flow</tt> in the energy balance. The default component uses a constant coefficient of heat transfer to model convective heat transfer between segment boundary (<tt>thermalPort</tt>) and the bulk flow. The <tt>heatTransfer</tt> model is replaceable and can be exchanged with any model extended from <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer\">PartialPipeHeatTransfer</a>. 
</html>",
      revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));

equation
//Pressure difference in momentum balance
dp=dp_stat;
  state=medium.state;

  //Pressure drop and gravity
if from_dp and not WallFriction.dp_is_zero then
    m_flow[1] = WallFriction.massFlowRate_dp(
      dp_stat[1] - height_ab/n*system.g*d[1],
      d_a,
      d[1],
      eta_a,
      eta[1],
      length/n,
      d_h,
      roughness,
      dp_small);
    for i in 2:n loop
      m_flow[i] = WallFriction.massFlowRate_dp(
        dp_stat[i] - height_ab/n*system.g*(d[i-1] + d[i])/2,
        d[i - 1],
        d[i],
        eta[i - 1],
        eta[i],
        length/n,
        d_h,
        roughness,
        dp_small);
    end for;
else
     dp_stat[1] = WallFriction.pressureLoss_m_flow(
        m_flow[1],
        d_a,
        d[1],
        eta_a,
        eta[1],
        length/n,
        d_h,
        roughness,
        m_flow_small) + height_ab/n*system.g*d[1];
    for i in 2:n loop
      dp_stat[i] = WallFriction.pressureLoss_m_flow(
        m_flow[i],
        d[i - 1],
        d[i],
        eta[i - 1],
        eta[i],
        length/n,
        d_h,
        roughness,
        m_flow_small) + height_ab/n*system.g*(d[i - 1] + d[i])/2;
    end for;
end if;
connect(thermalPort, heatTransfer.thermalPort) 
  annotation (Line(points={{0,54},{0,14}}, color={191,0,0}));
end DistributedPipeSb;

model DistributedPipeSa "Distributed pipe model"

  extends Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow(
    Qs_flow=heatTransfer.Q_flow,
    Ws_flow=zeros(n),
    ms_flow=zeros(n),
    msXi_flow=zeros(n, Medium.nXi),
    Vi=ones(n)*V/n,
    final modelStructure=Modelica_Fluid.Types.ModelStructure.av_b,
    redeclare final Modelica_Fluid.Interfaces.FluidStatePort_a port_a(redeclare
          package Medium = 
        Medium, m_flow(min=if allowFlowReversal and not static then -Modelica.Constants.inf else  0)),
    redeclare final Modelica_Fluid.Interfaces.FluidPort_b port_b(redeclare
          package Medium = 
        Medium, m_flow(max=if allowFlowReversal and not static then +Modelica.Constants.inf else 0)));
  parameter SI.Area area_h=perimeter*length "Heat transfer area" 
                                                             annotation(Dialog(tab="General", group="Heat transfer"));
  final parameter SI.Volume V=area*length;
  parameter SI.Length height_ab=0.0 "Height(port_b) - Height(port_a)" 
                                                                   annotation(Dialog(tab="General", group="Geometry"),Evaluate=true);

//Pressure Drop
 replaceable package WallFriction = 
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent
                                                                                constrainedby
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.PartialWallFriction
      "Characteristic of wall friction" 
                                       annotation(Dialog(tab="General", group="Pressure loss"),choicesAllMatching=true);
  parameter Boolean use_d_nominal=false
      "= true, if d_nominal is used, otherwise computed from medium"                              annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);
  parameter SI.Diameter d_h=4*area/perimeter "Hydraulic diameter" 
                                                                 annotation(Dialog(tab="General", group="Pressure loss"));
  parameter SI.Length roughness(min=0) = 2.5e-5
      "Absolute roughness of pipe (default = smooth steel pipe)" 
    annotation(Dialog(tab="General", group="Pressure loss",enable=WallFriction.use_roughness));
  parameter Boolean use_eta_nominal=false
      "= true, if eta_nominal is used, otherwise computed from medium"                          annotation(Dialog(tab="General", group="Pressure loss"),Evaluate=true);
  parameter SI.DynamicViscosity eta_nominal=Medium.dynamicViscosity(Medium.setState_pTX(Medium.p_default, Medium.T_default, Medium.X_default))
      "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
                                                                          annotation(Dialog(tab="General", group="Pressure loss",enable=use_nominal));
  parameter Boolean show_Re=false
      "= true, if Reynolds number is included for plotting" 
   annotation (Evaluate=true, Dialog(tab="Advanced", group="Pressure loss"));
  parameter Boolean from_dp=true
      " = true, use m_flow = f(dp), otherwise dp = f(m_flow)" 
  annotation (Evaluate=true, Dialog(tab="Advanced", group="Pressure loss"));
  parameter SI.AbsolutePressure dp_small=1
      "Turbulent flow if |dp| >= dp_small (only used if WallFriction=QuadraticTurbulent)"
  annotation(Dialog(tab="Advanced",group="Pressure loss", enable=from_dp and WallFriction.use_dp_small));

  parameter SI.MassFlowRate m_flow_small=0.01
      "Turbulent flow if |m_flow| >= m_flow_small (only used if WallFriction=QuadraticTurbulent)"
  annotation(Dialog(tab="Advanced",group="Pressure loss", enable=not from_dp and WallFriction.use_m_flow_small));

  SI.ReynoldsNumber[n+1] Re=Modelica_Fluid.Utilities.ReynoldsNumber_m_flow(
    m_flow,
    (cat(1, {eta_a}, eta) + cat(1, eta, {eta_b}))*0.5,
    diameter) if                                                                                               show_Re
      "Reynolds number of pipe flow";
  inner Medium.ThermodynamicState[n] state;
  replaceable Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_constAlpha
      heatTransfer(
    redeclare final package Medium = Medium,
    final n=n,
    final d_h=d_h,
    final A_h=area_h,
    final A_cross=area,
    final L=length,
    T=medium.T) constrainedby
      Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer(
    redeclare package Medium = Medium,
    n=n,
    d_h=d_h,
    A_h=area_h,
    A_cross=area,
    L=length,
    T=medium.T) "Convective heat transfer" 
            annotation (Dialog(tab="General", group="Heat transfer"),editButton=true,choicesAllMatching,
      Placement(transformation(extent={{-20,-20},{20,20}}, rotation=0)));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort
      "Thermal port" 
annotation (Placement(transformation(extent={{-10,44},{10,64}}, rotation=0)));
  outer Modelica_Fluid.System system "System properties";

  protected
  SI.Pressure[n] dp_stat;
   SI.DynamicViscosity eta_a=if not WallFriction.use_eta then 1.e-10 else (if 
      use_eta_nominal then eta_nominal else (if use_approxPortProperties then eta[1] else Medium.dynamicViscosity(Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)))));
  SI.DynamicViscosity eta_b=if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then eta_nominal else (if use_approxPortProperties then eta[n] else Medium.dynamicViscosity(Medium.setState_phX(port_b.p, inStream(port_b.h_outflow), inStream(port_b.Xi_outflow)))));
  SI.DynamicViscosity[n] eta=if not WallFriction.use_eta then fill(1.e-10, n) else (if use_eta_nominal then fill(eta_nominal, n) else 
      Medium.dynamicViscosity(medium.state));

annotation (
  Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={
          Rectangle(
            extent={{-100,44},{100,40}},
            lineColor={0,0,0},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-100,-40},{100,-44}},
            lineColor={0,0,0},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
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
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-150,-102},{150,-142}},
            lineColor={0,0,255},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255},
            textString="%name"),
          Line(
            points={{40,-75},{-50,-75}},
            color={0,128,255},
            smooth=Smooth.None),
          Polygon(
            points={{30,-60},{70,-75},{30,-90},{30,-60}},
            lineColor={0,128,255},
            smooth=Smooth.None,
            fillColor={0,128,255},
            fillPattern=FillPattern.Solid)}),
  Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}},
          grid={1,1}),
          graphics),
  Documentation(info="<html>
Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow\">PartialDistributedFlow</a>. Source terms in the mass balances are set to zero. The total volume is a parameter. The additional component <tt>heatTransfer</tt> specifies the source term <tt>Qs_flow</tt> in the energy balance. The default component uses a constant coefficient of heat transfer to model convective heat transfer between segment boundary (<tt>thermalPort</tt>) and the bulk flow. The <tt>heatTransfer</tt> model is replaceable and can be exchanged with any model extended from <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer\">PartialPipeHeatTransfer</a>. 
</html>",
      revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));

equation
//Pressure difference in momentum balance
dp=dp_stat;
  state=medium.state;

  //Pressure drop and gravity
if from_dp and not WallFriction.dp_is_zero then
    for i in 1:n-1 loop
      m_flow[i+1] = WallFriction.massFlowRate_dp(
        dp_stat[i] - height_ab/n*system.g*(d[i] + d[i+1])/2,
        d[i],
        d[i+1],
        eta[i],
        eta[i+1],
        length/n,
        d_h,
        roughness,
        dp_small);
    end for;
    m_flow[n+1] = WallFriction.massFlowRate_dp(
      dp_stat[n] - height_ab/n*system.g*d[n],
      d[n],
      d_b,
      eta[n],
      eta_b,
      length/n,
      d_h,
      roughness,
      dp_small);

else
    for i in 1:n-1 loop
      dp_stat[i] = WallFriction.pressureLoss_m_flow(
        m_flow[i+1],
        d[i],
        d[i+1],
        eta[i],
        eta[i+1],
        length/n,
        d_h,
        roughness,
        m_flow_small) + height_ab/n*system.g*(d[i] + d[i+1])/2;
    end for;
    dp_stat[n] = WallFriction.pressureLoss_m_flow(
      m_flow[n+1],
      d[n],
      d_b,
      eta[n],
      eta_b,
      length/n,
      d_h,
      roughness,
      m_flow_small) + height_ab/n*system.g*d[n];
end if;
connect(thermalPort, heatTransfer.thermalPort) 
  annotation (Line(points={{0,54},{0,14}}, color={191,0,0}));
end DistributedPipeSa;

  package BaseClasses
    extends Modelica_Fluid.Icons.BaseClassLibrary;

  partial model PartialDistributedFlowLumpedPressure
      "Base class for 1D fluid flow with the number of momentum balances reduced to 2"
    import Modelica_Fluid.Types;
    import Modelica.Constants.*;
    outer Modelica_Fluid.System system "System properties";
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
        "Fluid medium model" 
        annotation (choicesAllMatching=true);

  //Discretization
    parameter Integer n(min=1)=1 "Number of pipe segments";
    final parameter Integer np=2 "Number of momentum balances"                            annotation(Dialog(tab="Advanced"),Evaluate=true);
    final parameter Integer nl=integer(n/2)+1
        "Number of control volume that contains single state"                 annotation(Evaluate=true);

  //Advanced model options
    parameter Modelica_Fluid.Types.FlowDirection flowDirection=
        system.flowDirection
        "Unidirectional (port_a -> port_b) or bidirectional flow" 
       annotation(Dialog(tab="Advanced", group="Mass and energy balances"));
    parameter Boolean static=false
        "= true, static balances, no mass or energy is stored" 
                                  annotation(Dialog(tab="Advanced", group="Mass and energy balances"),Evaluate=true);

    parameter Boolean use_approxPortProperties=false
        "=true, port properties for pressure drop correlation are taken from neighboring control volume"
                                                                                                       annotation(Dialog(tab="Advanced", group="Momentum balances"),Evaluate=true);

  //Initialization
      parameter Types.Init initType=Types.
          Init.NoInit "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Boolean use_T_start=true
        "Use T_start if true, otherwise h_start" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p_a_start=Medium.p_default
        "Port a pressure start value" 
      annotation(Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p_b_start=p_a_start
        "Port b pressure start value" 
      annotation(Dialog(tab = "Initialization"));
      final parameter Medium.AbsolutePressure[n] p_start=if n > 1 then linspace(
          p_a_start - (p_a_start - p_b_start)/(2*n),
          p_b_start + (p_a_start - p_b_start)/(2*n),
          n) else {(p_a_start + p_b_start)/2} "Pressure start values";
      parameter Medium.Temperature T_start=if use_T_start then Medium.T_default else 
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
      parameter Medium.MassFlowRate mflow_start
        "Start value for mass flow rate"                                       annotation(Evaluate=true, Dialog(tab = "Initialization"));
      final parameter SI.Pressure dp_start=p_a_start - p_b_start;

  //Geometry parameters
    parameter Boolean isCircular=true
        "= true if cross sectional area is circular" 
      annotation (Evaluate, Dialog(tab="General", group="Geometry"));
    parameter SI.Length length "Length"   annotation(Dialog(tab="General", group="Geometry"));
    parameter SI.Diameter diameter "Diameter of circular pipe"      annotation(Dialog(group="Geometry", enable=isCircular));
    parameter SI.Length perimeter=if isCircular then Modelica.Constants.pi*diameter else 0
        "Inner perimeter"                                                                                     annotation(Dialog(tab="General", group="Geometry", enable=not isCircular));
    parameter SI.Area area=if isCircular then Modelica.Constants.pi*diameter*diameter/4 else 0
        "Inner cross section area"          annotation(Dialog(tab="General", group="Geometry", enable=not isCircular));
    input SI.Volume[n] Vi
        "Discretized volume, to be determined in inheriting class ";

  //Total quantities
    SI.Energy[n] U "Internal energy of fluid";
    SI.Mass[n] m "Mass of fluid";
    SI.Mass[n,Medium.nXi] mXi "Masses of independent components in the fluid";

  //Flow quantities
    inner Medium.MassFlowRate[n + 1] m_flow(each min=if allowFlowReversal then -inf else 
                0, each start=mflow_start, each fixed=false)
        "Mass flow rates of fluid across segment boundaries";
    SI.Velocity[n+1] v "velocity at volume boundaries (for display purposes)";
    Medium.MassFlowRate[n + 1,Medium.nXi] mXi_flow
        "Independent mass flow rates across segment boundaries";
    Medium.EnthalpyFlowRate[n + 1] H_flow
        "Enthalpy flow rates of fluid across segment boundaries";
   parameter Boolean use_d_nominal=false
        "= true, if d_nominal is used, otherwise computed from medium"                              annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);
   parameter SI.Density d_nominal=Medium.density_pTX(Medium.p_default, Medium.T_default, Medium.X_default)
        "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
                                                               annotation(Dialog(tab="Advanced", group="Momentum balance",enable=use_nominal));

    SI.Pressure[np] dp(start=dp0) "Pressure difference across staggered grid";

  //Fluid ports (the ports are replaceable, in order to be able to use different icons)
    replaceable Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare package
          Medium = 
          Medium, m_flow(min=if allowFlowReversal and not static then -inf else 0))
        "Fluid inlet port" 
                         annotation (Placement(transformation(extent={{-110,-10},
                {-90,10}}, rotation=0)));
    replaceable Modelica_Fluid.Interfaces.FluidPort_b port_b(redeclare package
          Medium = 
          Medium, m_flow(max=if allowFlowReversal and not static then +inf else 0))
        "Fluid outlet port" 
                          annotation (Placement(transformation(extent={{90,-10},
                {110,10}}, rotation=0)));
    Medium.BaseProperties[n] medium(
      each preferredMediumStates=if static then false else true,
      p(start=p_start),
      each h(start=h_start),
      each T(start=T_start),
      each Xi(start=X_start[1:Medium.nXi]));

     annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}),
                         graphics),
                          Icon(coordinateSystem(preserveAspectRatio=false,
              extent={{-100,-100},{100,100}}), graphics={Rectangle(
              extent={{-100,40},{100,-40}},
              lineColor={0,0,0},
              fillPattern=FillPattern.HorizontalCylinder,
              fillColor={0,127,255})}),
        Documentation(info="<html>
<p>The model <b>PartialDistributedFlow</b> is used as a base class for pipe flows with one-dimensional spatial discretization according to the finite volume method. The flow path is divided into <tt><b>n</b></tt> segments.</p>
 
<p><b>Mass and energy balances</b></p>
<p>One total mass and one energy balance is formed across each segment. If the medium contains more than one component, substance mass balances are added. Changes in potential and kinetic energy are neglected in the energy balance. The following source (or sink) terms are used in the balances and must be specified in extending models to complete this partial class:</p>
<ul>
<li>Energy balance: <tt><b>Qs_flow</b></tt>, e.g. convective or latent heat flow rate across segment boundary and <tt><b>Ws_flow</b></tt>, e.g. mechanical power</li>
<li>Total mass balance: <tt><b>ms_flow</b></tt>, e.g. condensing mass flow of negligible volume such as water in moist air</li>
<li>Substance mass balance: <tt><b>msXi_flow</b></tt>, as above</li>
</ul>
If the flag <tt>static</tt> is <b>true</b> then no mass or energy is stored in the component and the mass and energy balances are reduced to a quasi steady-state formulation. It should be noted that dynamic balances are required if flow reversal should be allowed.
<p>In addition the volume vector <tt><b>Vi</b></tt>, which specifies the volume of each segment and the pressure drop (or rise) in each segment <tt><b>dp</b></tt> must be provided in the extending class.
 
<p><b>Momentum balance</b></p>
<p>The momentum balance is always static, i.e. no dynamic momentum term is used. The momentum balances are formed across the segment boundaries (staggered grid). For this model only two momentum balances are formed on each side of a single pressure state (roughly half way along the flowpath). This assumes a constant pressure level for all medium models in the pipe. The total pressure drop (or rise) is split to be located on each end of the component. Connecting two pipes results in an algebraic pressure at the ports. Specifying a good start value for the port pressure is essential in order to solve large systems. The term <tt>dp</tt> is unspecified in this partial class. When extending from this model it may contain
<ul>
<li>pressure drop due to friction and other dissipative losses</li>
<li>changes in pressure resulting from significant variation of flow velocity along the flow path (with the assumption of a constant cross sectional area it must result from fluid density changes, such as in two-phase flow)</li>
<li>gravity effects for non-horizontal pipes</li>
</ul>
At least one relationship between pressure difference and massflow rate (dp=dp(m_flow)) per momentum balance is required in the extending class for this term to receive a fully determined model.
 
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

    //Source terms, have to be set in inheriting class
    protected
    parameter Boolean allowFlowReversal=
      flowDirection == Modelica_Fluid.Types.FlowDirection.Bidirectional
        "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
      annotation(Evaluate=true, Hide=true);
    input Medium.MassFlowRate[n] ms_flow "Mass flow rate, source or sink";
    input Medium.MassFlowRate[n,Medium.nXi] msXi_flow
        "Independent mass flow rates, source or sink";
    input SI.HeatFlowRate[n] Qs_flow
        "Heat and/or enthalpy flow rate, perpendicular to flow direction";
    input SI.Power[n] Ws_flow "Mechanical power, p*der(V) etc.";

    final parameter SI.Pressure[np] dp0=fill(dp_start/np,np);
    SI.Density[n] d=if use_d_nominal then ones(n)*d_nominal else medium.d;
    SI.Density d_a=if use_d_nominal then d_nominal else (if use_approxPortProperties then d[1] else Medium.density_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)));
    SI.Density d_b=if use_d_nominal then d_nominal else (if use_approxPortProperties then d[n] else Medium.density_phX(port_b.p, inStream(port_b.h_outflow), inStream(port_b.Xi_outflow)));

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
      v[i] = m_flow[i]/(medium[i - 1].d + medium[i].d)*2/area;
    end for;
    H_flow[1] = semiLinear(port_a.m_flow, inStream(port_a.h_outflow), medium[1].h);
    H_flow[n + 1] = -semiLinear(port_b.m_flow, inStream(port_b.h_outflow), medium[n].h);
    mXi_flow[1, :] = semiLinear(port_a.m_flow, inStream(port_a.Xi_outflow), medium[1].Xi);
    mXi_flow[n + 1, :] = -semiLinear(port_b.m_flow, inStream(port_b.Xi_outflow), medium[n].Xi);
    v[1] = m_flow[1]/(d_a + d[1])*2/area;
    v[n + 1] = m_flow[n + 1]/(d[n] + d_b)*2/area;

    // Total quantities
    for i in 1:n loop
      m[i] =Vi[i]*medium[i].d;
      mXi[i, :] = m[i]*medium[i].Xi;
      U[i] = m[i]*medium[i].u;
    end for;

    //Mass and energy balances
    if static then
    //steady state mass and energy balances, no numerical states, no flow reversal possible
      for i in 1:n loop
        0 = m_flow[i] - m_flow[i + 1] + ms_flow[i];
        zeros(Medium.nXi) = mXi_flow[i, :] - mXi_flow[i + 1, :] + msXi_flow[i, :];
        0 = H_flow[i] - H_flow[i + 1] + Qs_flow[i] + Ws_flow[i];
      end for;
    else
     //dynamic mass and energy balances, n "thermal" states, 1 pressure state
      for i in 1:n loop
        der(m[i]) = m_flow[i] - m_flow[i + 1] + ms_flow[i];
        der(mXi[i, :]) = mXi_flow[i, :] - mXi_flow[i + 1, :] + msXi_flow[i, :];
        der(U[i]) = H_flow[i] - H_flow[i + 1] + Qs_flow[i] + Ws_flow[i];
      end for;
      end if;
      for i in 1:n loop
        assert((allowFlowReversal and not static) or (m_flow[i] >= 0),
          "Flow reversal not allowed in distributed pipe");
      end for;

  //Momentum Balance, dp contains contributions from acceleration, gravitational and friction effects
  //two momentum balances, one on each side of pressure state, dp must be supplied in extending class
      dp[1] = port_a.p - medium[nl].p;
      dp[2] = medium[nl].p - port_b.p;
      if n == 2 then
        medium[2].p = medium[1].p;
      elseif n > 2 then
        medium[1:nl - 1].p = ones(nl - 1)*medium[nl].p;
        medium[nl + 1:n].p = ones(n - nl)*medium[nl].p;
      end if;

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
        if not Medium.singleState then
          der(medium[nl].p) = 0;
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
        medium[nl].p=p_start[nl];
        end if;
     elseif initType == Types.Init.SteadyStateHydraulic then
      //Steady state initialization for hydraulic states (p)
        if use_T_start then
          medium.T = ones(n)*T_start;
        else
          medium.h = ones(n)*h_start;
        end if;
        if not Medium.singleState then
        der(medium[nl].p) = 0;
        end if;
      else
        assert(false, "Unsupported initialization option");
      end if;
    end if;
  end PartialDistributedFlowLumpedPressure;

  partial model PartialDistributedFlow
    import Modelica_Fluid.Types;
    outer Modelica_Fluid.System system "System properties";
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
        "Fluid medium model" 
        annotation (choicesAllMatching=true);

  //Discretization
    parameter Integer n(min=1)=1 "Number of pipe segments";

  //Advanced model options
    parameter Modelica_Fluid.Types.FlowDirection flowDirection=
        system.flowDirection
        "Unidirectional (port_a -> port_b) or bidirectional flow" 
       annotation(Dialog(tab="Advanced", group="Mass and energy balances"));
    parameter Boolean static=false
        "= true, static balances, no mass or energy is stored" 
                                  annotation(Dialog(tab="Advanced", group="Mass and energy balances"),Evaluate=true);

    parameter Types.ModelStructure modelStructure=Types.ModelStructure.a_v_b
        "Determines whether flow model between volume and port is present"                                                                           annotation(Dialog(tab="Advanced", group="Mass and energy balances"),Evaluate=true);

    parameter Boolean use_approxPortProperties=false
        "=true, port properties for pressure drop correlation are taken from neighboring control volume"
                                                                                                       annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);

  //Initialization
      parameter Types.Init initType=Types.
          Init.NoInit "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Boolean use_T_start=true
        "Use T_start if true, otherwise h_start" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p_a_start=Medium.p_default
        "Start value of pressure at port a" 
      annotation(Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p_b_start=p_a_start
        "Start value of pressure at port b" 
      annotation(Dialog(tab = "Initialization"));
      final parameter Medium.AbsolutePressure[n] p_start=if n > 1 then linspace(
          p_a_start - (p_a_start - p_b_start)/(2*n),
          p_b_start + (p_a_start - p_b_start)/(2*n),
          n) else {(p_a_start + p_b_start)/2} "Start value of pressure";
      parameter Medium.Temperature T_start=if use_T_start then Medium.T_default else 
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
      parameter Medium.MassFlowRate mflow_start
        "Start value for mass flow rate" 
       annotation(Evaluate=true, Dialog(tab = "Initialization"));
      final parameter SI.Pressure dp_start=p_a_start - p_b_start;

  //Geometry
    parameter Boolean isCircular=true
        "= true if cross sectional area is circular" 
      annotation (Evaluate, Dialog(tab="General", group="Geometry"));
    parameter SI.Length length "Length"   annotation(Dialog(tab="General", group="Geometry"));
    parameter SI.Diameter diameter "Diameter of circular pipe"      annotation(Dialog(group="Geometry", enable=isCircular));
    parameter SI.Length perimeter=if isCircular then Modelica.Constants.pi*diameter else 0
        "Inner perimeter"                                                                                     annotation(Dialog(tab="General", group="Geometry", enable=not isCircular));
    parameter SI.Area area=if isCircular then Modelica.Constants.pi*diameter*diameter/4 else 0
        "Inner cross section area"          annotation(Dialog(tab="General", group="Geometry", enable=not isCircular));
    input SI.Volume[n] Vi "Discretized volume, determine in inheriting class ";

  //Total quantities
    SI.Energy[n] U "Internal energy of fluid";
    SI.Mass[n] m "Fluid mass";
    SI.Mass[n,Medium.nXi] mXi "Substance mass";

  //Flow quantities
    inner Medium.MassFlowRate[n + 1] m_flow(each min=if allowFlowReversal then -Modelica.Constants.inf else 
                0, each start=mflow_start, each fixed=false)
        "Mass flow rates of fluid across segment boundaries";
    SI.Velocity[n+1] v "Velocity at volume boundaries (not used in balances)";
    Medium.MassFlowRate[n + 1,Medium.nXi] mXi_flow
        "Independent mass flow rates across segment boundaries";
    Medium.EnthalpyFlowRate[n + 1] H_flow
        "Enthalpy flow rates of fluid across segment boundaries";
   parameter Boolean use_d_nominal=false
        "= true, if d_nominal is used, otherwise computed from medium"                              annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);
   parameter SI.Density d_nominal = Medium.density_pTX(Medium.p_default, Medium.T_default, Medium.X_default)
        "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
                                                               annotation(Dialog(tab="Advanced", group="Momentum balance",enable=use_d_nominal));
   SI.Pressure[np] dp(start=dp0) "pressure difference across staggered grid";

  //Fluid ports (the ports are replaceable, in order to be able to use different icons)
    replaceable Modelica_Fluid.Interfaces.FluidPort port_a "Fluid inlet port" 
                         annotation (Placement(transformation(extent={{-110,-10},
                {-90,10}}, rotation=0)));
    replaceable Modelica_Fluid.Interfaces.FluidPort port_b "Fluid outlet port" 
                          annotation (Placement(transformation(extent={{90,-10},
                {110,10}}, rotation=0)));
    Medium.BaseProperties[n] medium(
      each preferredMediumStates=if static then false else true,
      p(start=p_start),
      each h(start=h_start),
      each T(start=T_start),
      each Xi(start=X_start[1:Medium.nXi]));

     annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}),
                         graphics),
                          Icon(coordinateSystem(preserveAspectRatio=true,
              extent={{-100,-100},{100,100}}), graphics={Rectangle(
              extent={{-100,40},{100,-40}},
              lineColor={0,0,0},
              fillPattern=FillPattern.HorizontalCylinder,
              fillColor={0,127,255})}),
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
    //Source terms, have to be set in inheriting class (to zero if not used)
    protected
    parameter Boolean allowFlowReversal=
      flowDirection == Modelica_Fluid.Types.FlowDirection.Bidirectional
        "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
      annotation(Evaluate=true, Hide=true);
    input Medium.MassFlowRate[n] ms_flow "Mass flow rate, source or sink";
    input Medium.MassFlowRate[n,Medium.nXi] msXi_flow
        "Independent mass flow rates, source or sink";
    input SI.HeatFlowRate[n] Qs_flow "Heat flow rate, source or sink";
    input SI.Power[n] Ws_flow "Mechanical power, p*der(V) etc.";
    final parameter Integer np=if modelStructure==Types.ModelStructure.avb then n-1 else if (modelStructure==Types.ModelStructure.a_vb or modelStructure==Types.ModelStructure.av_b) then n else n+1;
    final parameter SI.Pressure[np] dp0=fill(dp_start,np)
        "pressure difference start values";
    SI.Density[n] d=if use_d_nominal then ones(n)*d_nominal else medium.d;
    SI.Density d_a=if use_d_nominal then d_nominal else (if use_approxPortProperties then d[1] else Medium.density_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)));
    SI.Density d_b=if use_d_nominal then d_nominal else (if use_approxPortProperties then d[n] else Medium.density_phX(port_b.p, inStream(port_b.h_outflow), inStream(port_b.Xi_outflow)));
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
      v[i] = m_flow[i]/(medium[i - 1].d + medium[i].d)*2/area;
    end for;
    H_flow[1] = semiLinear(port_a.m_flow, inStream(port_a.h_outflow), medium[1].h);
    H_flow[n + 1] = -semiLinear(port_b.m_flow, inStream(port_b.h_outflow), medium[n].h);
    mXi_flow[1, :] = semiLinear(port_a.m_flow, inStream(port_a.Xi_outflow), medium[1].Xi);
    mXi_flow[n + 1, :] = -semiLinear(port_b.m_flow, inStream(port_b.Xi_outflow), medium[n].Xi);
    v[1] = m_flow[1]/(d_a + d[1])*2/area;
    v[n + 1] = m_flow[n + 1]/(d[n] + d_b)*2/area;

    // Total quantities
    for i in 1:n loop
      m[i] =Vi[i]*medium[i].d;
      mXi[i, :] = m[i]*medium[i].Xi;
      U[i] = m[i]*medium[i].u;
    end for;

    //Mass and energy balances
    if static then
    //steady state mass and energy balances, no numerical states, no flow reversal possible
      for i in 1:n loop
        0 = m_flow[i] - m_flow[i + 1] + ms_flow[i];
        zeros(Medium.nXi) = mXi_flow[i, :] - mXi_flow[i + 1, :] + msXi_flow[i, :];
        0 = H_flow[i] - H_flow[i + 1] + Qs_flow[i];
      end for;
     else
    //dynamic mass and energy balances, n "thermal" states, n pressure states (if not singleState_hydraulic)
      for i in 1:n loop
        der(m[i]) = m_flow[i] - m_flow[i + 1] + ms_flow[i];
        der(mXi[i, :]) = mXi_flow[i, :] - mXi_flow[i + 1, :] + msXi_flow[i, :];
        der(U[i]) = H_flow[i] - H_flow[i + 1] + Qs_flow[i];
      end for;
    end if;
    for i in 1:n loop
      assert((allowFlowReversal and not static) or (m_flow[i] >= 0), "Flow reversal not allowed in distributed pipe");
    end for;

  //Simplified Momentum Balance, dp contains contributions from gravitational and friction effects
  //must be specified in extending class
  if modelStructure==Types.ModelStructure.av_b then
      port_a.p = medium[1].p;
      for i in 1:n-1 loop
        dp[i]=medium[i].p-medium[i+1].p;
      end for;
      dp[n]=medium[n].p-port_b.p;
  elseif modelStructure==Types.ModelStructure.a_vb then
      dp[1]=port_a.p-medium[1].p;
      for i in 2:n loop
        dp[i]=medium[i-1].p-medium[i].p;
      end for;
      port_b.p=medium[n].p;
  elseif modelStructure==Types.ModelStructure.avb then
      port_a.p=medium[1].p;
      for i in 1:n-1 loop
        dp[i]=medium[i].p-medium[i+1].p;
      end for;
      port_b.p=medium[n].p;
  else //ModelStructure.a_v_b
      dp[1]=port_a.p-medium[1].p;
      for i in 2:n loop
        dp[i]=medium[i-1].p-medium[i].p;
      end for;
      dp[n+1]=medium[n].p-port_b.p;
  end if;

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
  end PartialDistributedFlow;

    package CharacteristicNumbers
      function ReynoldsNumber
        input SI.MassFlowRate m_flow "Mass flow rate";
        input SI.Length d_ch "Characteristic length (hyd. diam. in pipes)";
        input SI.Area A "Cross sectional area";
        input SI.DynamicViscosity eta "Dynamic viscosity";
        output SI.ReynoldsNumber Re "Reynolds number";
        annotation (Documentation(info="Calculate Re-Number; Re = mdot*Dhyd/A/eta"),
             Icon(graphics));
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
      replaceable package Medium=Modelica.Media.Interfaces.PartialMedium annotation(Dialog(tab="No input", enable=false));
      parameter Integer n(min=1)=1 "Number of pipe segments" annotation(Dialog(tab="No input", enable=false));
      SI.HeatFlowRate[n] Q_flow "Heat flow rates";
      parameter SI.Area A_h "Total heat transfer area" annotation(Dialog(tab="No input", enable=false));
      parameter SI.Length d_h "Hydraulic diameter" annotation(Dialog(tab="No input", enable=false));
      parameter SI.Area A_cross "Cross flow area" annotation(Dialog(tab="No input", enable=false));
      parameter SI.Length L "Total pipe length" annotation(Dialog(tab="No input", enable=false));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort
          "Thermal port" 
        annotation (Placement(transformation(extent={{-20,60},{20,80}},
                rotation=0)));
      input SI.Temperature[n] T;
    equation

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
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
</html>"));
    end PartialPipeHeatTransfer;

    partial model PipeHT_Nu
        "base class for pipe heat transfer correlation in terms of Nusselt numberheat transfer in a circular pipe for laminar and turbulent one-phase flow"
      extends PartialPipeHeatTransfer;
      parameter SI.CoefficientOfHeatTransfer alpha0=100;
      outer input Medium.ThermodynamicState[n] state;
      outer input SI.MassFlowRate[n+1] m_flow;
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
      Re = CharacteristicNumbers.ReynoldsNumber(m_flow[1:n], fill(d_h,n), fill(A_cross, n), eta);
      Nu = CharacteristicNumbers.NusseltNumber(alpha, fill(d_h,n), lambda);
      thermalPort.Q_flow=Q_flow;
      for i in 1:n loop
        thermalPort[i].Q_flow=alpha[i]*A_h/n*(thermalPort[i].T-T[i]);
      end for;
        annotation (Documentation(info="<html>
Base class for heat transfer models that are expressed in terms of the Nusselt number and which can be used in distributed pipe models.
</html>"));
    end PipeHT_Nu;

    model PipeHT_constAlpha
      extends PartialPipeHeatTransfer;
      parameter SI.CoefficientOfHeatTransfer alpha0=200;
      annotation(Documentation(info="<html>
Simple heat transfer correlation with constant heat transfer coefficient, used as default component in <a distributed pipe models.
</html>"));
    equation
      for i in 1:n loop
        thermalPort[i].Q_flow=alpha0*A_h/n*(thermalPort[i].T-T[i]);
      end for;
      thermalPort.Q_flow=Q_flow;
    end PipeHT_constAlpha;
    annotation (Documentation(info="<html>
Heat transfer correlations for pipe models
</html>"));

    model PipeHT_LamTurb_local
        "laminar and turbulent forced convection in pipes, local coeff."
      extends PipeHT_Nu;
      protected
      Real[n] Nu_turb "Nusselt number for turbulent flow";
      Real[n] Nu_lam "Nusselt number for laminar flow";
      Real Nu_1;
      Real[n] Nu_2;
      Real[n] Xi;
      parameter SI.Length dx=L/n;
    equation
      Nu_1=3.66;
      for i in 1:n loop
       Nu_turb[i]=smooth(0,(Xi[i]/8)*abs(Re[i])*Pr[i]/(1+12.7*(Xi[i]/8)^0.5*(Pr[i]^(2/3)-1))*(1+1/3*(d_h/dx/(if m_flow[i]>=0 then (i-0.5) else (n-i+0.5)))^(2/3)));
       Xi[i]=(1.8*Modelica.Math.log10(max(1e-10,Re[i]))-1.5)^(-2);
       Nu_lam[i]=(Nu_1^3+0.7^3+(Nu_2[i]-0.7)^3)^(1/3);
       Nu_2[i]=smooth(0,1.077*(abs(Re[i])*Pr[i]*d_h/dx/(if m_flow[i]>=0 then (i-0.5) else (n-i+0.5)))^(1/3));
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
<li>The correlation holds for non-circular pipes only in the turbulent region. Use d_h=4*A/P as characteristic length.</li>
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
