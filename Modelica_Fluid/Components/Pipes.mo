package Pipes 
  
model LumpedPipe "Short pipe with one volume, wall friction and gravity effect" 
    import SI = Modelica.SIunits;
    
  extends Modelica_Fluid.Interfaces.Records.PartialInitializationParameters;
  replaceable package WallFriction = 
    Modelica_Fluid.SubClasses.PressureLosses.WallFriction.QuadraticTurbulent 
    extends Modelica_Fluid.Interfaces.PressureLosses.PartialWallFriction 
      "Characteristic of wall friction" 
                                       annotation(choicesAllMatching=true);
    
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
  parameter SI.DynamicViscosity eta_nominal = 0.01 
      "Nominal dynamic viscosity (for wall friction computation)" 
                                                                annotation(Dialog(enable=use_nominal));
  parameter SI.Density d_nominal = 0.01 
      "Nominal density (for wall friction computation)" 
                                                      annotation(Dialog(enable=use_nominal));
  parameter Types.FlowDirection.Temp flowDirection=
                  Types.FlowDirection.Unidirectional 
      "Unidirectional (port_a -> port_b) or bidirectional flow component" 
     annotation(Dialog(tab="Advanced"));
  parameter SI.AbsolutePressure dp_small = 1 
      "Turbulent flow if |dp| >= dp_small (only used if WallFriction=QuadraticTurbulent)"
    annotation(Dialog(tab="Advanced", enable=WallFriction.use_dp_small));
    
    Modelica_Fluid.Interfaces.Ports.FluidPort_a port_a(
                                  redeclare package Medium = Medium) 
      "Fluid connector a (positive design flow direction is from port_a to port_b)"
      annotation (extent=[-110,-10; -90,10]);
    Modelica_Fluid.Interfaces.Ports.FluidPort_b port_b(
                                  redeclare package Medium = Medium) 
      "Fluid connector b (positive design flow direction is from port_a to port_b)"
      annotation (extent=[110,-10; 90,10]);
    
  annotation (defaultComponentName="pipe",Icon(
      Rectangle(extent=[-100,44; 100,-44],   style(
          color=0,
          gradient=2,
          fillColor=8)),
      Rectangle(extent=[-100,40; 100,-40],   style(
          color=69,
          gradient=2,
          fillColor=69)),
      Text(
        extent=[-145,-40; 155,-90],
        string="%name",
        style(gradient=2, fillColor=69)),
      Ellipse(extent=[-11,10; 9,-10],   style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0}))),       Documentation(info="<html>
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
    Diagram,
    Coordsys(grid=[1,1], scale=0));
  PressureLosses.WallFrictionAndGravity frictionAndGravity1(
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
    show_Re=false)     annotation (extent=[-60,-10; -40,10]);
  Modelica_Fluid.SubClasses.ControlVolumes.PortVolume volume(
    redeclare package Medium = Medium,
    V=Modelica.Constants.pi*(diameter/2)^2*length,
    initType=initType,
    p_start=p_start,
    use_T_start=use_T_start,
    T_start=T_start,
    h_start=h_start,
    X_start=X_start) 
    annotation (extent=[-10,-10; 10,10]);
  PressureLosses.WallFrictionAndGravity frictionAndGravity2(
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
    show_Re=false)     annotation (extent=[40,-10; 60,10]);
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a thermalPort 
    annotation (extent=[-10,44; 10,64]);
equation 
  connect(frictionAndGravity1.port_a, port_a) 
    annotation (points=[-60,0; -100,0], style(color=69, rgbcolor={0,127,255}));
  connect(frictionAndGravity1.port_b, volume.port) 
    annotation (points=[-40,0; 0,0], style(color=69, rgbcolor={0,127,255}));
  connect(frictionAndGravity2.port_a, volume.port) 
    annotation (points=[40,0; 0,0], style(color=69, rgbcolor={0,127,255}));
  connect(frictionAndGravity2.port_b, port_b) 
    annotation (points=[60,0; 100,0], style(color=69, rgbcolor={0,127,255}));
  connect(volume.thermalPort, thermalPort) 
    annotation (points=[0,10; 0,54], style(color=42, rgbcolor={191,0,0}));
end LumpedPipe;
  
  model DistributedPipe_thermal "Distributed pipe model" 
    
    extends Modelica_Fluid.Interfaces.ControlVolumes.PartialDistributedFlow(
      Qs_flow=heatTransfer.Q_flow,
      ms_flow=zeros(n),
      msXi_flow=zeros(n, Medium.nXi),
      final singleState_thermal=false);
    parameter SI.Area area_h=perimeter*length "Heat transfer area" 
                                                               annotation(Dialog(tab="General", group="Heat transfer"));
    final parameter SI.Volume V=area*length;
    parameter SI.Length height_ab=0.0 "Height(port_b) - Height(port_a)" 
                                                                     annotation(Evaluate=true);
    
  //Pressure Drop
    parameter Boolean kineticTerm=false 
      " = true, include kinetic term in momentum balance"                                                                                 annotation(Evaluate=true);
    replaceable package WallFriction = 
        Modelica_Fluid.SubClasses.PressureLosses.WallFriction.QuadraticTurbulent
                                                                   extends 
      Modelica_Fluid.Interfaces.PressureLosses.PartialWallFriction 
      "Characteristic of wall friction"  annotation(Dialog(tab="General", group="Pressure loss"),choicesAllMatching=true);
    parameter Boolean use_d_nominal=false 
      "= true, if d_nominal is used, otherwise computed from medium"                                annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);
    parameter SI.Density d_nominal=0.01 
      "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
                                                               annotation(Dialog(tab="Advanced", group="Momentum balance",enable=use_nominal));
    parameter SI.Diameter d_h=4*area/perimeter "Hydraulic diameter" 
                                                                   annotation(Dialog(tab="General", group="Pressure loss"));
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
    SI.ReynoldsNumber[n] Re=Modelica_Fluid.Utilities.ReynoldsNumber_m_flow(
      m_flow,
      (eta_a + eta_b)/2,
      diameter) if                                                                                               show_Re 
      "Reynolds number of pipe flow";
    SI.Pressure[np] dp_stat;
    inner Medium.ThermodynamicState[n] state=medium.state 
      "may be used in heat transfer correlation";
    replaceable Modelica_Fluid.SubClasses.HeatTransfer.PipeHT_constAlpha 
      heatTransfer(
      redeclare final package Medium = Medium,
      final n=n,
      final d_h=d_h,
      final A_h=area_h,
      final A_cross=area,
      T=medium.T) extends 
      Modelica_Fluid.Interfaces.HeatTransfer.PartialPipeHeatTransfer(
      redeclare final package Medium = Medium,
      final n=n,
      final d_h=d_h,
      final A_h=area_h,
      final A_cross=area,
      T=medium.T) "Convective heat transfer" 
              annotation (Dialog(tab="General", group="Heat transfer"),choicesAllMatching, extent=[-20,-20;
      20,20]);
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort 
      "Thermal port" 
  annotation (extent=[-10,44; 10,64]);
    outer Components.Ambient ambient "Ambient conditions";
    
  protected 
    Real DI_flow[np] "Momentum flow (if kineticTerm=true)";
    SI.DynamicViscosity eta_a=if not WallFriction.use_eta then 1.e-10 else (if 
        use_eta_nominal then eta_nominal else eta[1]);
    SI.DynamicViscosity eta_b=if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then eta_nominal else eta[n]);
    SI.DynamicViscosity[n] eta=if not WallFriction.use_eta then fill(1.e-10, n) else (if use_eta_nominal then fill(eta_nominal, n) else 
        Medium.dynamicViscosity(medium));
    
  annotation (
    Icon(
      Rectangle(extent=[-100,44; 100,40], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95})),
      Rectangle(extent=[-100,-40; 100,-44], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95})),
      Ellipse(extent=[-72,10; -52,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Ellipse(extent=[-30,10; -10,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Ellipse(extent=[10,10; 30,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Ellipse(extent=[50,10; 70,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Text(
        extent=[-143,-42; 157,-92],
        string="%name",
        style(gradient=2, fillColor=69))),
    Diagram,
    Documentation(info="<html>
Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.BaseClasses.Pipes.PartialDistributedFlow\">PartialDistributedFlow</a>. Source terms in the mass balances are set to zero, and the option of reducing the thermal states to one is removed. All other model options remain. The additional component <tt>heatTransfer</tt> specifies the source term <tt>Qs_flow</tt> in the energy balance. The default component uses a constant coefficient of heat transfer to model convective heat transfer between segment boundary (<tt>thermalPort</tt>) and the bulk flow. The <tt>heatTransfer</tt> model is replaceable and can be exchanged with any model extended from <a href=\"Modelica:Modelica_Fluid.BaseClasses.Pipes.HeatTransfer.PartialPipeHeatTransfer\">PartialPipeHeatTransfer</a>. <b>DistributedPipe_thermal</b> is mainly designed for <b>thermal applications</b>, such as heat exchangers, where the transients of the internal energy may play an important role.
</html>",
        revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  equation 
    Vi = ones(n)*V/n;
    
  //Pressure difference in momentum balance
    dp = DI_flow + dp_stat;
    
    if kineticTerm and not singleState_hydraulic then
      for i in 2:n loop
        DI_flow[i] = Modelica_Fluid.Utilities.regSquare(m_flow[i], delta=max(
          0.0001, 0.01*mflow_start))/area*(1/d[i - 1] - 1/d[i]);
      end for;
      DI_flow[1] = Modelica_Fluid.Utilities.regSquare(m_flow[1], delta=max(0.0001,
        0.01*mflow_start))/area*(1/d_a - 1/d[1]);
      DI_flow[np] = Modelica_Fluid.Utilities.regSquare(m_flow[np], delta=max(
        0.0001, 0.01*mflow_start))/area*(1/d[n] - 1/d_b);
    else
      DI_flow = fill(0, np);
    end if;
    
    //Pressure drop and gravity
  if from_dp and not WallFriction.dp_is_zero then
    if singleState_hydraulic then
      m_flow[1] = WallFriction.massFlowRate_dp(
        dp[1] - ((integer(n/2) + 1)*2 - 1)/(2*n)*height_ab*ambient.g*(d_a + d_b)
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
        dp[2] - (2*n - (integer(n/2) + 1)*2 + 1)/(2*n)*height_ab*ambient.g*(d_a
           + d_b)/2,
        d_a,
        d_b,
        eta_a,
        eta_b,
        (2*n - (integer(n/2) + 1)*2 + 1)/(2*n)*length,
        d_h,
        roughness,
        dp_small);
    else
      m_flow[1] = WallFriction.massFlowRate_dp(
        dp[1] - height_ab/2/n*ambient.g*(d_a + d[1])/2,
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
          dp[i] - height_ab/n*ambient.g*(d[i - 1] + d[i])/2,
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
        dp[np] - height_ab/n/2*ambient.g*(d[n] + d_b)/2,
        d[n],
        d_b,
        eta[n],
        eta_b,
        length/n/2,
        d_h,
        roughness,
        dp_small);
    end if;
  else
    if singleState_hydraulic then
      dp[1] = (WallFriction.pressureLoss_m_flow(
        m_flow[1],
        d_a,
        d_b,
        eta_a,
        eta_b,
        ((integer(n/2) + 1)*2 - 1)/(2*n)*length,
        d_h,
        roughness,
        m_flow_small) + ((integer(n/2) + 1)*2 - 1)/(2*n)*height_ab*ambient.g*(
        d_a + d_b)/2);
      dp[2] = (WallFriction.pressureLoss_m_flow(
        m_flow[1],
        d_a,
        d_b,
        eta_a,
        eta_b,
        (2*n - (integer(n/2) + 1)*2 + 1)/(2*n)*length,
        d_h,
        roughness,
        m_flow_small) + (2*n - (integer(n/2) + 1)*2 + 1)/(2*n)*height_ab*
        ambient.g*(d_a + d_b)/2);
    else
      dp[1] = WallFriction.pressureLoss_m_flow(
        m_flow[1],
        d_a,
        d[1],
        eta_a,
        eta[1],
        length/n/2,
        d_h,
        roughness,
        m_flow_small) + height_ab/2/n*ambient.g*(d_a + d[1])/2;
      for i in 2:n loop
        dp[i] = WallFriction.pressureLoss_m_flow(
          m_flow[i],
          d[i - 1],
          d[i],
          eta[i - 1],
          eta[i],
          length/n,
          d_h,
          roughness,
          m_flow_small) + height_ab/n*ambient.g*(d[i - 1] + d[i])/2;
      end for;
      dp[np] = WallFriction.pressureLoss_m_flow(
        m_flow[np],
        d[n],
        d_b,
        eta[n],
        eta_b,
        length/n/2,
        d_h,
        roughness,
        m_flow_small) + height_ab/n/2*ambient.g*(d[n] + d_b)/2;
    end if;
  end if;
  connect(thermalPort, heatTransfer.thermalPort) 
    annotation (points=[0,54; 0,14], style(color=42, rgbcolor={191,0,0}));
  end DistributedPipe_thermal;
  
 model DistributedPipe_hydraulic "Distributed pipe model" 
    
   extends Modelica_Fluid.Interfaces.ControlVolumes.PartialDistributedFlow(
     Qs_flow=zeros(n),
     ms_flow=zeros(n),
     msXi_flow=zeros(n, Medium.nXi),
     final singleState_hydraulic=false);
    
    //Volume size
   final parameter SI.Volume V=area*length;
    
    //Pressure Drop
  parameter Boolean kineticTerm=false 
      " = true, include kinetic term in momentum balance"            annotation(Evaluate=true);
  replaceable package WallFriction = 
      Modelica_Fluid.SubClasses.PressureLosses.WallFriction.QuadraticTurbulent 
    extends Modelica_Fluid.Interfaces.PressureLosses.PartialWallFriction 
      "Characteristic of wall friction" 
                                       annotation(Dialog(tab="General", group="Pressure loss"),choicesAllMatching=true);
  parameter SI.Diameter d_h=4*area/perimeter "Hydraulic diameter" annotation(Dialog(tab="General", group="Pressure loss"));
   parameter SI.Length height_ab=0.0 "Height(port_b) - Height(port_a)" 
                                                                     annotation(Evaluate=true);
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
   SI.ReynoldsNumber[n] Re=Modelica_Fluid.Utilities.ReynoldsNumber_m_flow(
       m_flow,
       (eta_a + eta_b)/2,
       diameter) if                                                                                              show_Re 
      "Reynolds number of pipe flow";
   SI.Pressure[np] dp_stat;
   outer Components.Ambient ambient "Ambient conditions";
    
  protected 
   Real DI_flow[np];
   SI.DynamicViscosity eta_a=if not WallFriction.use_eta then 1.e-10 else (if 
       use_eta_nominal then eta_nominal else eta[1]);
    //approximation
   SI.DynamicViscosity eta_b=if not WallFriction.use_eta then 1.e-10 else (if 
       use_eta_nominal then eta_nominal else eta[n]);                                                                      //approximation
   SI.DynamicViscosity[n] eta=if not WallFriction.use_eta then fill(1.e-10, n) else (if use_eta_nominal then fill(eta_nominal, n) else 
       Medium.dynamicViscosity(medium));
    
    annotation (
      Icon(
        Rectangle(extent=[-100,44; 100,40], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95})),
                         Rectangle(extent=[-100,-40; 100,-44], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95})),
      Ellipse(extent=[-72,10; -52,-10], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0})),
      Ellipse(extent=[-30,10; -10,-10], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0})),
      Ellipse(extent=[10,10; 30,-10],   style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0})),
      Ellipse(extent=[50,10; 70,-10],   style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0})),
      Text(
        extent=[-143,-42; 157,-92],
        string="%name",
        style(gradient=2, fillColor=69))),
                          Diagram,
    Documentation(info="<html>
Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.BaseClasses.Pipes.PartialDistributedFlow\">PartialDistributedFlow</a>. Source terms in mass and energy balances are set to zero, and the option of reducing the pressure states to one is removed. All other model options remain. <b>DistributedPipe_hydraulic</b> is mainly designed for <b>hydraulic applications</b>, such a long insulated pipes, where the pressure transients may play an important role but not thermal behaviour. If no significant effect on medium properties along the flow path is expected consider using the <a href=\"Modelica:Modelica_Fluid.Components.Pipes.LumpedPipe\">LumpedPipe</a> model, which uses one set of dynamic medium states only or <a href=\"Modelica:Modelica_Fluid.Components.PressureLosses.WallFrictionAndGravity\">PressureLosses.WallFrictionAndGravity</a>, without storage of energy and mass. For thermal applications, especially heat transfer across the pipe wall, use instead <a href=\"Modelica:Modelica_Fluid.Components.Pipes.DistributedPipe_thermal\">DistributedPipe_thermal</a>.
</html>",
        revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
    
 equation 
   Vi = ones(n)*V/n;
    
 //Pressure difference in momentum balance
   dp = DI_flow + dp_stat;
    
   //Momentum flow terms
   if kineticTerm then
     for i in 2:n loop
       DI_flow[i] = Modelica_Fluid.Utilities.regSquare(m_flow[i], delta=max(
         0.0001, 0.01*mflow_start))/area*(1/d[i - 1] - 1/d[i]);
     end for;
     DI_flow[1] = Modelica_Fluid.Utilities.regSquare(m_flow[1], delta=max(0.0001,
       0.01*mflow_start))/area*(1/d_a - 1/d[1]);
     DI_flow[np] = Modelica_Fluid.Utilities.regSquare(m_flow[np], delta=max(
       0.0001, 0.01*mflow_start))/area*(1/d[n] - 1/d_b);
   else
     DI_flow = fill(0, np);
   end if;
    
    //Friction loss and gravity
   if from_dp and not WallFriction.dp_is_zero then
     m_flow[1] = WallFriction.massFlowRate_dp(
       dp_stat[1] - height_ab/2/n*ambient.g*(d_a + d[1])/2, d_a, d[1], eta_a, eta[1], length/n/2, d_h,
       roughness, dp_small);
     for i in 2:n loop
       m_flow[i] = WallFriction.massFlowRate_dp(
         dp_stat[i] - height_ab/n*ambient.g*(d[i - 1] + d[i])/2, d[i - 1], d[i], eta[i - 1], eta[i],
         length/n, d_h, roughness, dp_small);
     end for;
     m_flow[n + 1] = WallFriction.massFlowRate_dp(
       dp_stat[np] - height_ab/n/2*ambient.g*(d[n] + d_b)/2, d[n], d_b, eta[n], eta_b, length/n/2,
       d_h, roughness, dp_small);
   else
     dp_stat[1] = WallFriction.pressureLoss_m_flow(
       m_flow[1], d_a, d[1], eta_a, eta[1], length/n/2, d_h, roughness, m_flow_small) + height_ab/2/n*ambient.g*(d_a + d[1])/2;
     for i in 2:n loop
       dp_stat[i] = WallFriction.pressureLoss_m_flow(m_flow[i], d[i - 1], d[i], eta[i - 1], eta[i],
         length/n, d_h, roughness, m_flow_small) + height_ab/n*ambient.g*(d[i - 1] + d[i])/2;
     end for;
     dp_stat[np] = WallFriction.pressureLoss_m_flow(
       m_flow[np], d[n], d_b, eta[n], eta_b, length/n/2, d_h, roughness, m_flow_small) + height_ab/n/2*ambient.g*(d[n] + d_b)/2;
   end if;
 end DistributedPipe_hydraulic;
end Pipes;
