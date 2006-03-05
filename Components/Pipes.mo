package Pipes 
  
model LumpedPipe "Short pipe with one volume, wall friction and gravity effect" 
    import SI = Modelica.SIunits;
    
  extends BaseClasses.Common.PartialInitializationParameters;
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
    
  extends BaseClasses.Pipes.PartialDistributedFlow(
    Qs_flow=heatTransfer.Q_flow,
    ms_flow=zeros(n),
    msXi_flow=zeros(n, Medium.nXi),
    final singleState_thermal=false);
    
  parameter SI.Area area_h = perimeter*length "Heat transfer area" 
                                                                 annotation(Dialog(tab="General", group="Heat transfer"));
  inner Medium.ThermodynamicState[n] state = medium.state;
    
  replaceable BaseClasses.Pipes.HeatTransfer.PipeHT_constAlpha heatTransfer(
      redeclare final package Medium = Medium,
      final n=n,
      final d_h=d_h,
      final A_h=area_h,
      T=medium.T) extends 
      BaseClasses.Pipes.HeatTransfer.PartialPipeHeatTransfer(
      redeclare final package Medium = Medium,
      final n=n,
      final d_h=d_h,
      final A_h=area_h,
      T=medium.T) "Convective heat transfer" 
                annotation (Dialog(tab="General", group="Heat transfer"),choicesAllMatching, extent=[-20,-20;
        20,20]);
    
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort 
      "Thermal port" 
    annotation (extent=[-10,44; 10,64]);
  annotation (Icon(Rectangle(extent=[-100,44; 100,40], style(
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
Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.BaseClasses.Pipes.PartialDistributedFlow\">PartialDistributedFlow</a>. Source terms in the mass balances are set to zero, and the option of reducing the thermal states to one is removed. All other model options remain. The additional component <tt>heatTransfer</tt> specifies the source term <tt>Qs_flow</tt> in the energy balance. The default component uses a constant coefficient of heat transfer to model convective heat transfer between segment boundary (<tt>thermalPort</tt>) and the bulk flow. The <tt>heatTransfer</tt> model is replaceable and can be exchanged with any model extended from <a href=\"Modelica:Modelica_Fluid.BaseClasses.Pipes.HeatTransfer.PartialPipeHeatTransfer\">PartialPipeHeatTransfer</a>. <b>DistributedPipe_thermal</b> is mainly designed for <b>thermal applications</b>, such as heat exchangers, where the transients of the internal energy may play an important role.
</html>", revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  equation 
    
    connect(thermalPort, heatTransfer.thermalPort) 
      annotation (points=[0,54; 0,14], style(color=42, rgbcolor={191,0,0}));
  end DistributedPipe_thermal;

  model DistributedPipe_hydraulic "Distributed pipe model" 
    
  extends BaseClasses.Pipes.PartialDistributedFlow(
    Qs_flow=zeros(n),
    ms_flow=zeros(n),
    msXi_flow=zeros(n, Medium.nXi),
    final singleState_hydraulic=false);
    
  annotation (Icon(Rectangle(extent=[-100,44; 100,40], style(
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
Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.BaseClasses.Pipes.PartialDistributedFlow\">PartialDistributedFlow</a>. Source terms in mass and energy balances are set to zero, and the option of reducing the pressure states to one is removed. All other model options remain. <b>DistributedPipe_hydraulic</b> is mainly designed for <b>hydraulic applications</b>, such a long insulated pipes, where the transients of the pressure may play an important role but not thermal behaviour.
</html>", revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  end DistributedPipe_hydraulic;
end Pipes;
