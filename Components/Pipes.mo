package Pipes 
  model DistributedPipeFV_mitWand "Distributed pipe model with optional wall" 
    
  extends BaseClasses.Pipes.PartialFlow1D_FV(
    Qs_flow=heat.Q_flow,
    ms_flow=zeros(n),
    msXi_flow=zeros(n, Medium.nXi));
    
  parameter SI.Area area_h = P_inner*length "Heat transfer area" annotation(Dialog(tab="General", group="Heat transfer"));
    
  parameter Boolean use_wall=false 
      "= true, use wall component between fluid and thermalPort" 
                                                                annotation(Dialog(tab="General", group="Wall - optional"),Evaluate=true);
  parameter SI.Diameter d_outer "Outer diameter of circular pipe" annotation(Dialog(tab="General", group="Wall - optional", enable=(crossSectionType==1 and use_wall)));
  parameter SI.Length h_outer "Outer height of rectangular pipe"  annotation(Dialog(tab="General", group="Wall - optional", enable=(crossSectionType==2 and use_wall)));
  parameter SI.Length w_outer "Outer width of rectangular pipe"  annotation(Dialog(tab="General", group="Wall - optional", enable=(crossSectionType==2 and use_wall)));
  parameter SI.Length A_outer = if crossSectionType == 1 then Modelica.Constants.pi*d_outer*d_outer/4 else if crossSectionType == 2 then h_outer*w_outer else 1 
      "Outer cross section area" 
                               annotation(Dialog(tab="General", group="Wall - optional", enable=(use_wall and crossSectionType==3)));
  inner Medium.ThermodynamicState[n] state = medium.state;
    
  replaceable BaseClasses.Pipes.HeatTransfer.PipeHT_constAlpha heat(
    redeclare final package Medium = Medium,
    final n=n,
    final d_h=d_h,
    final A_h=area_h,
    T=medium.T) extends BaseClasses.Pipes.HeatTransfer.PartialPipeHeatTransfer(
    redeclare final package Medium = Medium,
    final n=n,
    final d_h=d_h,
    final A_h=area_h,
    T=medium.T) "Convective heat transfer" 
                annotation (Dialog(tab="General", group="Heat transfer"),choicesAllMatching, extent=[-20,-20;
        20,20]);
    
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort 
      "Thermal port" 
    annotation (extent=[-20,60; 20,80]);
  replaceable model Wall = 
        Components.Thermal.WallConstProps  extends 
      Components.Thermal.PartialPipeWall "Wall model"              annotation(choicesAllMatching, Dialog(enable=use_wall, tab="General", group="Wall - optional"));
  Wall wall(final n=n, final a_inner=A_inner, final a_outer=A_outer, final 
        length=length) if use_wall 
                           annotation (extent=[10,20; 50,60]);
  annotation (Icon(Rectangle(extent=[-100,60; 100,40], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95},
          fillPattern=7)), Rectangle(extent=[-100,-40; 100,-60], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95},
          fillPattern=7)),
      Text(
        extent=[-100,-60; 100,-100],
        string="%name",
        style(color=3, rgbcolor={0,0,255}))),
                            Diagram,
      Documentation(info="<html>
<p>
From Katrins email, Nov. 28, 2005:
</p>
 
<p>
extends Interfaces.1DFlow. Pressure drop and heat transfer are added in terms of replaceable components. The main problem here is to make all required variables and parameters available to the respective component (medium state, pipe geometry, Medium functions, empirical parameters). Only those shared by all future replaceable models (the simple one parameter model and the highly sophisticated (fictitious) two phase Nusselt correlation) can be set by modifiers (which is not straightforward in Dymola at the moment if a contsraining clause is used).  Those not required by all models as i.e. viscosity and conductivitiy must be computed inside the component from medium properties made available via inner and outer. I always try to avoid this as it it as bit like free climbing, but in this case I see no better solution.
</p>
 
<p>
Martin, I have not tested your latest pressure drop implementation with this model, but will do so as soon as possible. However, it is used in a completely different way, that means as an array of components, not as a  base class, in order to be able to handle distributed flow. I will check if another implementation would be more practical.
</p>
 
<p>
The pipe model contains a Boolean flag useWall which determines if a wall component is added. Unfortunately the icon does not represent the difference. In this way a heat exchanger can be created using two instances of the pipe model, one with a wall and one without. If interested in transients it could also make sense to include a wall in an insulated pipe. 
</p>
 
</html>"));
  equation 
    
    if use_wall then
  connect(wall.thermalPort_a, thermalPort) annotation (points=[30,50; 30,60; 0,
          60; 0,70],
      style(
      color=10,
      rgbcolor={95,95,95},
      fillColor=0,
      rgbfillColor={0,0,0},
      fillPattern=7));
  connect(wall.thermalPort_b, heat.thermalPort) 
                                              annotation (points=[30,30; 30,22;
          0,22; 0,14],
      style(
      color=10,
      rgbcolor={95,95,95},
      fillColor=0,
      rgbfillColor={0,0,0},
      fillPattern=7));
    else
  connect(thermalPort, heat.thermalPort) 
                                       annotation (points=[0,70; 0,14],
      style(
      color=10,
      rgbcolor={95,95,95},
      fillColor=0,
      rgbfillColor={0,0,0},
      fillPattern=7));
    end if;
  end DistributedPipeFV_mitWand;
  
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
is totally isolated (= no thermal flow from the fluid to the
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

  model DistributedPipeFV "Distributed pipe model" 
    
  extends BaseClasses.Pipes.PartialFlow1D_FV(
    Qs_flow=heatTransfer.Q_flow,
    ms_flow=zeros(n),
    msXi_flow=zeros(n, Medium.nXi));
    
  parameter SI.Area area_h = P_inner*length "Heat transfer area" annotation(Dialog(tab="General", group="Heat transfer"));
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
<p>
From Katrins email, Nov. 28, 2005:
</p>
 
<p>
extends Interfaces.1DFlow. Pressure drop and heat transfer are added in terms of replaceable components. The main problem here is to make all required variables and parameters available to the respective component (medium state, pipe geometry, Medium functions, empirical parameters). Only those shared by all future replaceable models (the simple one parameter model and the highly sophisticated (fictitious) two phase Nusselt correlation) can be set by modifiers (which is not straightforward in Dymola at the moment if a contsraining clause is used).  Those not required by all models as i.e. viscosity and conductivitiy must be computed inside the component from medium properties made available via inner and outer. I always try to avoid this as it it as bit like free climbing, but in this case I see no better solution.
</p>
 
<p>
Martin, I have not tested your latest pressure drop implementation with this model, but will do so as soon as possible. However, it is used in a completely different way, that means as an array of components, not as a  base class, in order to be able to handle distributed flow. I will check if another implementation would be more practical.
</p>
 
<p>
The pipe model contains a Boolean flag useWall which determines if a wall component is added. Unfortunately the icon does not represent the difference. In this way a heat exchanger can be created using two instances of the pipe model, one with a wall and one without. If interested in transients it could also make sense to include a wall in an insulated pipe. 
</p>
 
</html>"));
  equation 
    
    connect(thermalPort, heatTransfer.thermalPort)
      annotation (points=[0,54; 0,14], style(color=42, rgbcolor={191,0,0}));
  end DistributedPipeFV;
end Pipes;
