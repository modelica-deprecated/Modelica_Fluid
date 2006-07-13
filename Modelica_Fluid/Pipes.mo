package Pipes 
    extends Modelica_Fluid.Icons.VariantLibrary;
  
model LumpedPipe "Short pipe with one volume, wall friction and gravity effect" 
    import SI = Modelica.SIunits;
    
  extends Modelica_Fluid.Interfaces.PartialInitializationParameters;
  replaceable package WallFriction = 
    BaseClasses.PressureLosses.WallFriction.QuadraticTurbulent 
    extends BaseClasses.PressureLosses.WallFriction.PartialWallFriction 
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
    
    Interfaces.FluidPort_a port_a(redeclare package Medium = Medium) 
      "Fluid connector a (positive design flow direction is from port_a to port_b)"
      annotation (extent=[-110,-10; -90,10]);
    Interfaces.FluidPort_b port_b(redeclare package Medium = Medium) 
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
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity frictionAndGravity1(
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
      show_Re=false)   annotation (extent=[-60,-10; -40,10]);
  BaseClasses.Pipes.PortVolume volume(
    redeclare package Medium = Medium,
    V=Modelica.Constants.pi*(diameter/2)^2*length,
    initType=initType,
    p_start=p_start,
    use_T_start=use_T_start,
    T_start=T_start,
    h_start=h_start,
    X_start=X_start) 
    annotation (extent=[-10,-10; 10,10]);
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity frictionAndGravity2(
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
      show_Re=false)   annotation (extent=[40,-10; 60,10]);
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
    
    extends BaseClasses.Pipes.PartialDistributedFlow(
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
        BaseClasses.PressureLosses.WallFriction.QuadraticTurbulent extends 
      BaseClasses.PressureLosses.WallFriction.PartialWallFriction 
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
    replaceable BaseClasses.Pipes.HeatTransfer.PipeHT_constAlpha heatTransfer(
      redeclare final package Medium = Medium,
      final n=n,
      final d_h=d_h,
      final A_h=area_h,
      final A_cross=area,
      T=medium.T) extends 
      BaseClasses.Pipes.HeatTransfer.PartialPipeHeatTransfer(
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
    outer Modelica_Fluid.Ambient ambient "Ambient conditions";
    
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
        dp_stat[1] - ((integer(n/2) + 1)*2 - 1)/(2*n)*height_ab*ambient.g*(d_a + d_b)
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
        dp_stat[2] - (2*n - (integer(n/2) + 1)*2 + 1)/(2*n)*height_ab*ambient.g*(d_a
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
        dp_stat[1] - height_ab/2/n*ambient.g*(d_a + d[1])/2,
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
          dp_stat[i] - height_ab/n*ambient.g*(d[i - 1] + d[i])/2,
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
        dp_stat[np] - height_ab/n/2*ambient.g*(d[n] + d_b)/2,
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
      dp_stat[1] = (WallFriction.pressureLoss_m_flow(
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
      dp_stat[2] = (WallFriction.pressureLoss_m_flow(
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
      dp_stat[1] = WallFriction.pressureLoss_m_flow(
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
        dp_stat[i] = WallFriction.pressureLoss_m_flow(
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
      dp_stat[np] = WallFriction.pressureLoss_m_flow(
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
    
   extends BaseClasses.Pipes.PartialDistributedFlow(
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
      BaseClasses.PressureLosses.WallFriction.QuadraticTurbulent 
    extends BaseClasses.PressureLosses.WallFriction.PartialWallFriction 
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
   outer Modelica_Fluid.Ambient ambient "Ambient conditions";
    
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
  
  package BaseClasses 
  partial model PartialDistributedFlow 
      import Modelica_Fluid.Types;
      import Modelica.Constants.*;
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Fluid medium model" 
        annotation (choicesAllMatching=true);
      
  //Discretization
    parameter Integer n(min=1)=1 "Number of pipe segments";
    final parameter Integer np=if singleState_hydraulic then 2 else n + 1 
        "Number of momentum balances"                                                     annotation(Dialog(tab="Advanced"),Evaluate=true);
    final parameter Integer nl=integer(n/2)+1 
        "Number of control volume that contains single state"                 annotation(Evaluate=true);
      
  //Advanced model options
    parameter Boolean allowFlowReversal=true 
        "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
                                                                       annotation(Dialog(tab="Advanced", group="Mass and energy balances", enable=not static));
    parameter Boolean static=false 
        "= true, static balances, no mass or energy is stored" 
                                  annotation(Dialog(tab="Advanced", group="Mass and energy balances"),Evaluate=true);
    parameter Boolean singleState_hydraulic=false 
        " = true, lumped pressure drop, reduces number of pressure states to one"
                                                                                annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);
    parameter Boolean singleState_thermal=false 
        " = true, number of temperature or enthalpy states is reduced to one" 
                                                                             annotation(Evaluate=true,Dialog(tab="Advanced", group="Mass and energy balances", enable=not static));
      
  //Initialization
      parameter Types.Init.Temp initType=Types.
          Init.NoInit "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Boolean use_T_start=true 
        "Use T_start if true, otherwise h_start" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p_a_start=Medium.p_default 
        "Start value of pressure at port a" 
      annotation(Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p_b_start=Medium.p_default 
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
        "Start value for mass flow rate"                                       annotation(Evaluate=true, Dialog(tab = "Initialization"));
      final parameter SI.Pressure dp_start=p_a_start - p_b_start;
      
  //Geometry parameters
    parameter Modelica_Fluid.Types.CrossSectionTypes.Temp crossSectionType=
        Modelica_Fluid.Types.CrossSectionTypes.Circular 
        "Type of cross section of pipe" 
      annotation (Dialog(tab="General", group="Geometry"));
    parameter SI.Length length "Length"   annotation(Dialog(tab="General", group="Geometry"));
    parameter SI.Diameter diameter "Diameter of circular pipe"      annotation(Dialog(group="Geometry", enable=crossSectionType==1));
    parameter SI.Length height "Height of rectangular pipe"           annotation(Dialog(group="Geometry", enable=crossSectionType==2));
    parameter SI.Length width "Width of rectangular pipe"            annotation(Dialog(group="Geometry", enable=crossSectionType==2));
    parameter SI.Length perimeter=if crossSectionType == 1 then Modelica.Constants.pi*diameter else if crossSectionType == 2 then 2*height + 2*
        width else 1 "Inner perimeter"  annotation(Dialog(tab="General", group="Geometry", enable=crossSectionType==3));
    inner parameter SI.Area area=if crossSectionType == 1 then Modelica.Constants.pi*diameter*diameter/4 else if crossSectionType
         == 2 then height*width else 1 "Inner cross section area" 
                                            annotation(Dialog(tab="General", group="Geometry", enable=crossSectionType==3));
    SI.Volume[n] Vi "Discretized volume, determine in inheriting class ";
      
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
   parameter SI.Density d_nominal=0.01 
        "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
                                                               annotation(Dialog(tab="Advanced", group="Momentum balance",enable=use_nominal));
      
  //Source terms, have to be set in inheriting class
    Medium.MassFlowRate[n] ms_flow "Mass flow rate, source or sink";
    Medium.MassFlowRate[n,Medium.nXi] msXi_flow 
        "Independent mass flow rates, source or sink";
    SI.HeatFlowRate[n] Qs_flow "Heat flow rate, source or sink";
    SI.Pressure[np] dp(start=dp0) "pressure difference across staggered grid";
      
  //Fluid ports
    Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = 
          Medium, m_flow(min=if allowFlowReversal and not static then -inf else 0)) 
        "Fluid inlet port" 
                         annotation (extent=[-110,-10; -90,10]);
    Modelica_Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = 
          Medium, m_flow(max=if allowFlowReversal and not static then +inf else 0)) 
        "Fluid outlet port" 
                          annotation (extent=[90,-10; 110,10]);
    Medium.BaseProperties[n] medium(
      each preferredMediumStates=if static then false else true,
      p(start=p_start),
      each h(start=h_start),
      each T(start=T_start),
      each Xi(start=X_start[1:Medium.nXi]));
      
     annotation (Diagram, Icon(Rectangle(extent=[-100,40; 100,-40], style(
            color=69,
            gradient=2,
            fillColor=69))),
        Documentation(info="<html>
<p>The model <b>PartialDistributedFlow</b> is used as a base class for pipe flows with one-dimensional spatial discretization according to the finite volume method, such as <a href=\"Modelica:Modelica_Fluid.Components.Pipes.DistributedPipe_thermal\">DistributedPipe_thermal</a> and <a href=\"Modelica:Modelica_Fluid.Components.Pipes.DistributedPipe_hydraulic\">DistributedPipe_hydraulic</a>. The flow path is divided into <tt><b>n</b></tt> segments.</p>
 
<p><b>Mass and energy balances</b></p>
<p>One total mass and one energy balance is formed across each segment. If the medium contains more than one component, substance mass balances are added. Changes in potential and kinetic energy are neglected in the energy balance. The following source (or sink) terms are used in the balances and must be specified in extending models to complete this partial class:</p>
<ul>
<li>Energy balance: <tt><b>Qs_flow</b></tt>, e.g. convective or latent heat flow rate across segment boundary</li>
<li>Total mass balance: <tt><b>ms_flow</b></tt>, e.g. condensing mass flow of negligible volume such as water in air</li>
<li>Substance mass balance: <tt><b>msXi_flow</b></tt>, as above</li>
</ul>
<p>In addition the volume vector <tt><b>Vi</b></tt>, which specifies the volume of each segment and the pressure drop (or rise) in each segment <tt><b>dp</b></tt> must be provided in the extending class.
 
<p><b>Momentum balance</b></p>
<p>The momentum balance is always static, i.e. no dynamic momentum term is used. The momentum balances are formed across the segment boundaries (staggered grid). Half a momentum balance is used on each end of the flow model resulting in a total of n-1 full and 2 half momentum balances. Connecting two pipes therefore results in an algebraic pressure at the ports. Specifying a good start value for the port pressure is essential in order to solve large systems. The term <tt>dp</tt> is unspecified in this partial class. When extending from this model it may contain
<ul>
<li>pressure drop due to friction and other dissipative losses</li>
<li>changes in pressure resulting from significant variation of flow velocity along the flow path (with the assumption of a constant cross sectional area it must result from fluid density changes, such as in two-phase flow)</li>
<li>gravity effects for non-horizontal pipes</li>
</ul>
At least one relationship between pressure difference and massflow rate (dp=dp(m_flow)) is required in the extending class for this term to receive a fully determined model.
 
When connecting two components, e.g. two pipes, the momentum balance across the connection point reduces to</p> 
<pre>pipe1.port_b.p = pipe2.port_a.p</pre>
<p>This is only true if the flow velocity remains the same on each side of the connection. For any significant change in diameter (and if the resulting effects, such as change in kinetic energy, cannot be neglected) an adapter component should be used. This also allows taking into account friction losses with respect to the actual geometry of the connection point.</p>
 
<p><b>Reducing the number of numerical states</b></p>
<p>The default settings of this model result in <tt>n*(2 + nX - 1)</tt> numerical states with <tt>nX</tt> as the number of substances: <tt>n</tt> pressure states, <tt>n</tt> enthalpy or temperature states and <tt>n*nX-1</tt> mass fraction states. Depending on the simulation task the model efficiency may be increased if the number of numerical states is reduced. The following model options exist:</p>
<ul>
<li><tt><b>static</b></tt> - if true, no numerical states are present, static mass and energy balances</li>
<li><tt><b>singleState_hydraulic</b></tt> - if true, only one pressure state is present, just two momentum balances, algebraic constraints for remaining pressures</li>
<li><tt><b>singleState_thermal</b></tt> - if true, only one enthalpy or temperature state is present, just one dynamic energy balance, algebraic constraints for remaining properties</li>
</ul>
Selecting a fixed medium composition along the entire flow path is currently not possible.
<p>
 
</html>",   revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
    protected 
    final parameter SI.Pressure[np] dp0={if singleState_hydraulic then dp_start else dp_start/(if i
         > 1 and i < np then n else 2*n) for i in 1:np};
    SI.Density[n] d=if use_d_nominal then ones(n)*d_nominal else medium.d;
    /*SI.Density d_a=if use_d_nominal then d_nominal else Medium.density_ph(port_a.p, port_a.h);
  SI.Density d_b=if use_d_nominal then d_nominal else Medium.density_ph(port_b.p, port_b.h);
  Above equations currently produce nonlinear systems and numerical Jacobians*/
    SI.Density d_a=if use_d_nominal then d_nominal else d[1];//approximation
    SI.Density d_b=if use_d_nominal then d_nominal else d[n];//approximation
      
  equation 
    // Boundary conditions
    port_a.H_flow = semiLinear(port_a.m_flow, port_a.h, medium[1].h);
    port_b.H_flow = semiLinear(port_b.m_flow, port_b.h, medium[n].h);
    port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, medium[1].Xi);
    port_b.mXi_flow = semiLinear(port_b.m_flow, port_b.Xi, medium[n].Xi);
    port_a.m_flow = m_flow[1];
    port_b.m_flow = -m_flow[n + 1];
      
    // Distributed flow quantities
    for i in 2:n loop
      H_flow[i] = semiLinear(m_flow[i], medium[i - 1].h, medium[i].h);
      mXi_flow[i, :] = semiLinear(m_flow[i], medium[i - 1].Xi, medium[i].Xi);
      v[i] = m_flow[i]/(medium[i - 1].d + medium[i].d)*2/area;
    end for;
    H_flow[1] = port_a.H_flow;
    H_flow[n + 1] = -port_b.H_flow;
    mXi_flow[1, :] = port_a.mXi_flow;
    mXi_flow[n + 1, :] = -port_b.mXi_flow;
    v[1] = m_flow[1]/d_a/area;
    v[n + 1] = m_flow[n + 1]/d_b/area;
      
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
    elseif singleState_thermal then
    //dynamic mass balances, one dynamic energy balance, n pressure states (if not singleState_hydraulic), 1 "thermal" (h or T) state
      for i in 1:n loop
        der(m[i]) = m_flow[i] - m_flow[i + 1] + ms_flow[i];
        der(mXi[i, :]) = mXi_flow[i, :] - mXi_flow[i + 1, :] + msXi_flow[i, :];
      end for;
      der(U[nl]) = H_flow[nl] - H_flow[nl + 1] + Qs_flow[nl];
      for i in 1:nl - 1 loop
        medium[i].h = medium[nl].h;
      end for;
      for i in nl + 1:n loop
        medium[i].h = medium[nl].h;
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
      
  //Momentum Balance, dp contains contributions from acceleration, gravitational and friction effects
  if singleState_hydraulic then //two momentum balances, one on each side of pressure state
      dp[1] = port_a.p - medium[nl].p;
      dp[2] = medium[nl].p - port_b.p;
      if n == 2 then
        medium[2].p = medium[1].p;
      elseif n > 2 then
        medium[1:nl - 1].p = ones(nl - 1)*medium[nl].p;
        medium[nl + 1:n].p = ones(n - nl)*medium[nl].p;
      end if;
    else
      dp[1]=port_a.p-medium[1].p;
      for i in 2:n loop
        dp[i]=medium[i-1].p-medium[i].p;
      end for;
      dp[np]=medium[n].p-port_b.p;
    end if;
      
  initial equation 
    // Initial conditions
    if not static then
      if initType == Types.Init.NoInit then
      // no initial equations
      elseif initType == Types.Init.SteadyState then
      //steady state initialization
        if singleState_thermal then
          if use_T_start then
          der(medium[nl].T) = 0;
        else
          der(medium[nl].h) = 0;
        end if;
        else
        if use_T_start then
          der(medium.T) = zeros(n);
        else
          der(medium.h) = zeros(n);
        end if;
        end if;
        if not (singleState_hydraulic or Medium.singleState) then
          der(medium.p) = zeros(n);
        elseif singleState_hydraulic then
         der(medium[nl].p) = 0;
        end if;
        for i in 1:n loop
          der(medium[i].Xi) = zeros(Medium.nXi);
        end for;
      elseif initType == Types.Init.InitialValues then
      //Initialization with initial values
        if singleState_thermal then
          if use_T_start then
          medium[nl].T = T_start;
        else
          medium[nl].h = h_start;
        end if;
        else
        if use_T_start then
          medium.T = ones(n)*T_start;
        else
          medium.h = ones(n)*h_start;
        end if;
        end if;
        if not (singleState_hydraulic or Medium.singleState) then
           medium.p=p_start;
        elseif singleState_hydraulic then
         medium[nl].p=p_start[nl];
        end if;
      elseif initType == Types.Init.SteadyStateHydraulic then
      //Steady state initialization for hydraulic states (p)
        if use_T_start then
          medium.T = ones(n)*T_start;
        else
          medium.h = ones(n)*h_start;
        end if;
        if not (singleState_hydraulic or Medium.singleState) then
          der(medium.p) = zeros(n);
        elseif singleState_hydraulic then
          der(medium[nl].p) = 0;
        end if;
      else
        assert(false, "Unsupported initialization option");
      end if;
    end if;
  end PartialDistributedFlow;
    
  model PortVolume 
      "Fixed volume associated with a port by the finite volume method (used to build up physical components; fulfills mass and energy balance)" 
      import Modelica_Fluid.Types;
    extends Modelica_Fluid.Interfaces.PartialInitializationParameters;
      
    replaceable package Medium = 
        Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching = true);
      
    parameter SI.Volume V "Volume";
      
    Interfaces.FluidPort_a port(
      redeclare package Medium = Medium) "Fluid port" 
      annotation (extent=[-10, -10; 10, 10], rotation=0);
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a thermalPort 
        "Thermal port" 
      annotation (extent=[-10,90; 10,110]);
      
    Medium.BaseProperties medium(preferredMediumStates=true,
                p(start=p_start), T(start=T_start),
                h(start=h_start), Xi(start= X_start[1:Medium.nXi]));
    SI.Energy U "Internal energy of fluid";
    SI.Mass m "Mass of fluid";
    SI.Mass mXi[Medium.nXi] "Masses of independent components in the fluid";
    annotation (
     Icon(
        Ellipse(extent=[-100, 100; 100, -100], style(
            color=0,
            rgbcolor={0,0,0},
            gradient=3,
            fillColor=68,
            rgbfillColor={170,213,255})),
        Text(extent=[-150,-100; 150,-150], string="%name")),
                           Documentation(info="<html>
<p>
This component models the <b>volume</b> of <b>fixed size</b> that is
associated with the <b>fluid port</b> to which it is connected.
This means that all medium properties inside the volume, are identical
to the port medium properties. In particular, the specific enthalpy
inside the volume (= medium.h) is always identical to the specific enthalpy
in the port (port.h = medium.h). Usually, this model is used when
discretizing a component according to the finite volume method into
volumes in internal ports that only store energy and mass and into
transport elements that just transport energy, mass and momentum
between the internal ports without storing these quantities during the
transport. This splitting is only possible under certain assumptions.
</p>
</html>"),
      Diagram);
  equation 
    // medium properties set by port values
      port.p = medium.p;
      port.h = medium.h;
      port.Xi = medium.Xi;
      thermalPort.T = medium.T;
      
    // Total quantities
       m    = V*medium.d "Total Mass";
       mXi = m*medium.Xi "Independent component masses";
       U    = m*medium.u "Internal energy";
      
    // Mass and energy balance
       der(m)    = port.m_flow "Total mass balance";
       der(mXi)  = port.mXi_flow "Independent component mass balance";
       der(U)    = port.H_flow + thermalPort.Q_flow "Energy balance";
      
  initial equation 
    // Initial conditions
    if initType == Types.Init.NoInit then
      // no initial equations
    elseif initType == Types.Init.InitialValues then
      if not Medium.singleState then
         medium.p = p_start;
      end if;
      if use_T_start then
        medium.T = T_start;
      else
        medium.h = h_start;
      end if;
      medium.Xi = X_start[1:Medium.nXi];
    elseif initType == Types.Init.SteadyState then
      if not Medium.singleState then
         der(medium.p) = 0;
      end if;
      der(medium.h) = 0;
      der(medium.Xi) = zeros(Medium.nXi);
    elseif initType == Types.Init.SteadyStateHydraulic then
      if not Medium.singleState then
         der(medium.p) = 0;
      end if;
      if use_T_start then
        medium.T = T_start;
      else
        medium.h = h_start;
      end if;
      medium.Xi = X_start[1:Medium.nXi];
    else
      assert(false, "Unsupported initialization option");
    end if;
  end PortVolume;
    
  package HeatTransfer 
    partial model PartialPipeHeatTransfer 
      replaceable package Medium=Modelica.Media.Interfaces.PartialMedium annotation(Dialog(tab="No input", enable=false));
      parameter Integer n(min=1)=1 "Number of pipe segments" annotation(Dialog(tab="No input", enable=false));
      SI.HeatFlowRate[n] Q_flow "Heat flow rates";
      parameter SI.Area A_h "Total heat transfer area" annotation(Dialog(tab="No input", enable=false));
      parameter SI.Length d_h "Hydraulic diameter" annotation(Dialog(tab="No input", enable=false));
      parameter SI.Area A_cross "Cross flow area" annotation(Dialog(tab="No input", enable=false));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort 
          "Thermal port" 
        annotation (extent=[-20,60; 20,80]);
      SI.Temperature[n] T;
    equation 
        
      annotation (Icon(Ellipse(extent=[-60,64; 60,-56], style(
              color=42,
              rgbcolor={127,0,0},
              gradient=3,
              fillColor=1,
              rgbfillColor={232,0,0})), Text(
            extent=[-38,26; 40,-14],
            style(
              color=42,
              rgbcolor={127,0,0},
              gradient=3,
              fillColor=1,
              rgbfillColor={232,0,0},
              fillPattern=7),
            string="%name")), Documentation(info="<html>
Base class for heat transfer models that can be used in model <a href=\"Modelica:Modelica_Fluid.Components.Pipes.DistributedPipe_thermal\">DistributedPipe_thermal</a>.
</html>"));
    end PartialPipeHeatTransfer;
      
    model PipeHT_constAlpha 
      extends PartialPipeHeatTransfer;
      parameter SI.CoefficientOfHeatTransfer alpha0=200;
      annotation(structurallyIncomplete=true, Documentation(info="<html>
Simple heat transfer correlation with constant heat transfer coefficient, used as default component in <a href=\"Modelica:Modelica_Fluid.Components.Pipes.DistributedPipe_thermal\">DistributedPipe_thermal</a>.
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
  end HeatTransfer;
  end BaseClasses;
end Pipes;
