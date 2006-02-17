package Components "Basic components for fluid models" 
  
model FluidOptions 
    "Default options and environment settings for components of Modelica_Fluid" 
    import SI = Modelica.SIunits;
    import Modelica_Fluid.Types.FlowDirection;
    import Modelica_Fluid.Types.Init;
    import Modelica.SIunits.Conversions;
    
  parameter Init.Temp default_initOption = Init.NoInit 
      "Default initialization option" 
    annotation(Dialog(group="Defaults"));
  parameter FlowDirection.Temp default_flowDirection=FlowDirection.Bidirectional 
      "Default flow direction defined via Advanced.flowDirection" 
    annotation(Dialog(group="Defaults"));
  parameter Modelica.Media.Interfaces.PartialMedium.AbsolutePressure 
      default_p_ambient =                                                                101325 
      "Default ambient pressure" 
      annotation(Dialog(group="Defaults"));
  parameter Modelica.Media.Interfaces.PartialMedium.Temperature 
      default_T_ambient=
      Conversions.from_degC(20) "Default ambient temperature" 
      annotation(Dialog(group="Defaults"));
  parameter SI.Acceleration g=9.81 "Constant gravity acceleration" annotation(Dialog(group="Environment"));
    
  annotation (
    preferedView="info",
    defaultComponentName="fluidOptions",
    defaultComponentPrefixes="inner",
    missingInnerMessage="An inner \"fluidOptions\" component is not defined. A default 
fluidOptions component will be used. If this is not desired, 
drag Modelica_Fluid.Components.FluidOptions into the top level of your model.",
    Icon(
      Rectangle(extent=[-100,100; 100,-100], style(
          color=3,
          rgbcolor={0,0,255},
          fillColor=7,
          rgbfillColor={255,255,255})),
      Text(
        extent=[-160,160; 160,110],
        style(color=3, rgbcolor={0,0,255}),
        string="%name"),
      Line(points=[-86,-30; 82,-30], style(color=0, rgbcolor={0,0,0})),
      Line(points=[-82,-68; -52,-30], style(color=0, rgbcolor={0,0,0})),
      Line(points=[-48,-68; -18,-30], style(color=0, rgbcolor={0,0,0})),
      Line(points=[-14,-68; 16,-30], style(color=0, rgbcolor={0,0,0})),
      Line(points=[22,-68; 52,-30], style(color=0, rgbcolor={0,0,0})),
      Line(points=[74,84; 74,14], style(color=0, rgbcolor={0,0,0})),
      Polygon(points=[60,14; 88,14; 74,-18; 60,14], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Text(
        extent=[16,20; 60,-18],
        style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0},
          fillPattern=1),
        string="g"),
        Text(
          extent=[-90,82; 74,50],
          style(color=0, rgbcolor={0,0,0}),
          string="defaults")),
    Diagram,
    Documentation(info="<HTML>
<p>
This models defines <b>default options</b> (such as the default 
for \"allowFlowReversal\") for all components that are on the same 
or on a lower level as this component, as well as the constant 
gravity acceleration. Dragging this component in a model results
in the following declaration:
</p>
<pre>
   <b>inner</b> Modelica_Fluid.Components.FluidOptions fluidOptions;
</pre>
<p>
The parameters of this instance can be 
then accessed via a corresponding outer declaration:
</p>
<pre>
   <b>outer</b> Modelica_Fluid.Components.FluidOptions fluidOptions;
</pre>
<p>
Note, all parameters under group \"Defaults\" are used as 
default setting by the Modelica_Fluid components. They can
be individually redefined in the corresponding
component.
</p>
</HTML>
"));
    
    Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape dummyShape 
      "Just temporarily, to force Dymola to open an animation window (only then animation setup is available for diagram animation)"
      annotation (extent=[-60,20; -40,40]);
end FluidOptions;
  
  model MixingVolume 
    "Mixing volume with inlet and outlet ports (flow reversal is allowed)" 
    import Modelica_Fluid.Types.Init;
    import Modelica.Constants;
    import Modelica_Fluid.Types.FlowDirection;
    import Modelica_Fluid.Types.FlowDirectionWithGlobalDefault;
    extends Interfaces.PartialInitializationParameters;
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching = true);
    parameter SI.Volume V "Volume";
    parameter FlowDirectionWithGlobalDefault.Temp flowDirection=
              FlowDirectionWithGlobalDefault.UseGlobalFluidOption 
      "Unidirectional (port_a -> port_b) or bidirectional flow component" 
       annotation(Dialog(tab="Advanced"));
    Interfaces.FluidPort_a port_a(redeclare package Medium = Medium,
                       m_flow(min=if allowFlowReversal then -Constants.inf else 0)) 
      "Fluid inlet port" annotation (extent=[-112,-10; -92,10]);
    Interfaces.FluidPort_b port_b(redeclare package Medium = Medium,
                       m_flow(max=if allowFlowReversal then +Constants.inf else 0)) 
      "Fluid outlet port" annotation (extent=[90,-10; 110,10]);
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a thermalPort 
      "Thermal port" 
      annotation (extent=[-20,88; 20,108]);
    Medium.BaseProperties medium(preferredMediumStates=true,
                 p(start=p_start), h(start=h_start),
                 T(start=T_start), Xi(start=X_start[1:Medium.nXi]));
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
        Text(extent=[-144, 178; 146, 116], string="%name"),
        Text(
          extent=[-130, -108; 144, -150],
          style(color=0),
          string="V=%V")), Documentation(info="<html>
</html>"),
      Diagram);
  protected 
    outer Modelica_Fluid.Components.FluidOptions fluidOptions 
      "Global default options";
    parameter Boolean allowFlowReversal=
       flowDirection == FlowDirectionWithGlobalDefault.Bidirectional
       or flowDirection == FlowDirectionWithGlobalDefault.UseGlobalFluidOption
       and fluidOptions.default_flowDirection ==FlowDirection.Bidirectional 
      "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
       annotation(Evaluate=true, Hide=true);
  equation 
    // boundary conditions
    port_a.p = medium.p;
    port_b.p = medium.p;
    port_a.H_flow = semiLinear(port_a.m_flow, port_a.h, medium.h);
    port_b.H_flow = semiLinear(port_b.m_flow, port_b.h, medium.h);
    port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, medium.Xi);
    port_b.mXi_flow = semiLinear(port_b.m_flow, port_b.Xi, medium.Xi);
    thermalPort.T = medium.T;
    
    // Total quantities
    m    = V*medium.d;
    mXi = m*medium.Xi;
    U    = m*medium.u;
    
    // Mass and energy balance
    der(m)   = port_a.m_flow + port_b.m_flow;
    der(mXi) = port_a.mXi_flow + port_b.mXi_flow;
    der(U)   = port_a.H_flow + port_b.H_flow + thermalPort.Q_flow;
    
  initial equation 
    // Initial conditions
    if initOption == Init.NoInit then
      // no initial equations
    elseif initOption == Init.InitialValues then
      if not Medium.singleState then
        medium.p = p_start;
      end if;
      if use_T_start then
        medium.T = T_start;
      else
        medium.h = h_start;
      end if;
      medium.Xi = X_start[1:Medium.nXi];
    elseif initOption == Init.SteadyState then
      if not Medium.singleState then
         der(medium.p) = 0;
      end if;
      der(medium.h) = 0;
      der(medium.Xi) = zeros(Medium.nXi);
    elseif initOption == Init.SteadyStateHydraulic then
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
  end MixingVolume;
  
  model ShortPipe 
    "Short pipe with one volume, wall friction and gravity effect" 
    import SI = Modelica.SIunits;
    
    extends Interfaces.PartialInitializationParameters;
    extends Modelica_Fluid.Interfaces.PartialTwoPort;
    replaceable package WallFriction = 
      Modelica_Fluid.PressureLosses.Utilities.WallFriction.QuadraticTurbulent 
      extends 
      Modelica_Fluid.PressureLosses.Utilities.WallFriction.PartialWallFriction 
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
    parameter SI.DynamicViscosity eta_nominal = 0.01 
      "Nominal dynamic viscosity (for wall friction computation)" annotation(Dialog(enable=use_nominal));
    parameter SI.Density d_nominal = 0.01 
      "Nominal density (for wall friction computation)" annotation(Dialog(enable=use_nominal));
    parameter Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.Temp 
      flowDirection=
              Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.UseGlobalFluidOption 
      "Unidirectional (port_a -> port_b) or bidirectional flow component" 
       annotation(Dialog(tab="Advanced"));
    parameter SI.AbsolutePressure dp_small = 1 
      "Turbulent flow if |dp| >= dp_small (only used if WallFriction=QuadraticTurbulent)"
      annotation(Dialog(tab="Advanced", enable=WallFriction.use_dp_small));
    
    annotation (defaultComponentName="pipe",Icon(
        Rectangle(extent=[-100,60; 100,-60],   style(
            color=0,
            gradient=2,
            fillColor=8)),
        Rectangle(extent=[-100,34; 100,-36],   style(
            color=69,
            gradient=2,
            fillColor=69)),
        Text(
          extent=[-150,-60; 150,-110],
          string="%name",
          style(gradient=2, fillColor=69)),
        Ellipse(extent=[-15,15; 15,-15],  style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=0,
              rgbfillColor={0,0,0})),
        Line(points=[0,0; 0,70], style(color=42, rgbcolor={194,0,0}))),
                                             Documentation(info="<html>
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
    Utilities.PortVolume volume(
      redeclare package Medium = Medium,
      V=Modelica.Constants.pi*(diameter/2)^2*length,
      initOption=initOption,
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
      annotation (extent=[-10,60; 10,80]);
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
      annotation (points=[0,10; 0,70], style(color=42, rgbcolor={191,0,0}));
  end ShortPipe;
  
  extends Modelica.Icons.Library;
  
  annotation (preferedView="info",
              Documentation(info="<html>
<p>
This package will contain basic component models of the fluid library.
It is currently empty as all components are being evaluated together 
with examples first (see sub-package Examples). 
</p>
</html>"));
  
  model StaticHead 
    "Models the static head between two ports at different heights" 
    extends Modelica_Fluid.Interfaces.PartialTwoPortTransport;
    parameter SI.Height H_b_a "Height of port b over port a";
    parameter SI.Acceleration g = Modelica.Constants.g_n "Gravity acceleration"
                                                                                 annotation(Dialog(tab="Advanced"));
    Medium.Density d "Fluid density";
    annotation (Icon(
        Rectangle(extent=[-100,60; 100,-60],   style(
            color=0,
            gradient=2,
            fillColor=8)),
        Rectangle(extent=[-100,34; 100,-36],   style(
            color=69,
            gradient=2,
            fillColor=69)),
        Text(
          extent=[-150,140; 150,80],
          string="%name",
          style(gradient=2, fillColor=69))), Documentation(info="<html>
<p>
This model describes the static head due to the relative height between the two connectors. No mass, energy and momentum storage, and no pressure drop due to friction are considered.
</p>
</html>", revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Added to Modelica_Fluid</li>
</ul>
</html>"));
  equation 
   d = if dp > 0 then medium_a.d else medium_b.d;
   port_a.p = port_b.p + H_b_a*g*d;
  end StaticHead;
  
  model OpenTank "Tank with three inlet/outlet-arrays at variable heights" 
    import SI = Modelica.SIunits;
    import Modelica_Fluid.Types;
    replaceable package Medium = PackageMedium 
      extends Modelica.Media.Interfaces.PartialMedium "Medium in the component"
      annotation (choicesAllMatching=true);
    
    parameter SI.Height height "Height of tank";
    parameter SI.Area area "Area of tank";
    parameter SI.Volume V0=0 "Volume of the liquid when the level is zero";
    // parameter Real k = 4.9 "Heat transfer coefficient from tank to ambient";
    
    parameter Integer n_topPorts = 0 "Number of topPorts" 
      annotation(Dialog(group="topPorts (= pipes at top of tank; only flow into tank)"));
    parameter SI.Height top_heights[n_topPorts]=fill(height,n_topPorts) 
      "Heights of topPorts" 
      annotation(Dialog(group="topPorts (= pipes at top of tank; only flow into tank)",enable=n_topPorts > 0));
    
    parameter Integer n_bottomPorts = 0 "Number of bottomPorts" 
       annotation(Dialog(group="bottomPorts (= pipes at bottom of tank; in and out flow of tank)"));
    parameter SI.Height bottom_heights[n_bottomPorts]=fill(0.0,n_bottomPorts) 
      "Heights of bottomPorts" 
       annotation(Dialog(group="bottomPorts (= pipes at bottom of tank; in and out flow of tank)",enable=n_bottomPorts > 0));
    parameter SI.Diameter bottom_diameters[n_bottomPorts] 
      "Inner (hydraulic) diameters of bottomPorts" 
       annotation(Dialog(group="bottomPorts (= pipes at bottom of tank; in and out flow of tank)", enable=n_bottomPorts > 0));
    
    parameter Integer n_sidePorts = 0 "Number of sidePorts" 
       annotation(Dialog(group="sidePorts (= pipes at side of tank; in and out flow of tank)"));
    parameter SI.Height side_heights[n_sidePorts] "Heights of sidePorts" 
       annotation(Dialog(group="sidePorts (= pipes at side of tank; in and out flow of tank)",enable=n_sidePorts > 0));
    parameter SI.Area side_diameters[n_sidePorts] 
      "Inner (hydraulic) diameters of sidePorts" 
       annotation(Dialog(group="sidePorts (= pipes at side of tank; in and out flow of tank)",enable=n_sidePorts > 0));
    
    parameter Medium.AbsolutePressure p_ambient=fluidOptions.default_p_ambient 
      "Tank surface pressure" 
      annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));
    parameter Medium.Temperature T_ambient = fluidOptions.default_T_ambient 
      "Tank surface Temperature" 
      annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));
    
    parameter Types.InitWithGlobalDefault.Temp initOption=
              Types.InitWithGlobalDefault.InitialValues "Initialization option"
      annotation(Dialog(tab = "Ambient and Initialization", group = "Initialization"));
    parameter SI.Height level_start "Start value of tank level" 
      annotation(Dialog(tab="Ambient and Initialization", group = "Initialization"));
    parameter Medium.Temperature T_start=T_ambient "Start value of temperature"
      annotation(Dialog(tab = "Ambient and Initialization", group = "Initialization"));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Ambient and Initialization", group = "Initialization", enable=Medium.nXi > 0));
    
  /*
  parameter Real k_small(min=0) = 0 
    "Small regularization range if tank level is below bottom_heights[i] or side_heights[i]; k_small = 0 gives ideal switch"
              annotation(Evaluate=true, Dialog(tab="Advanced"));
*/
    
    Modelica_Fluid.Interfaces.FluidPort_ArrayIcon topPorts[n_topPorts](redeclare 
        package Medium = Medium, m_flow(each start=0), mXi_flow(each start=0)) 
      annotation (extent=[-30,100; 30,108]);
    Modelica_Fluid.Interfaces.FluidPort_ArrayIcon bottomPorts[n_bottomPorts](redeclare 
        package Medium = Medium, m_flow(each start=0), mXi_flow(each start=0)) 
      annotation (extent=[-30,-108; 30,-100],   rotation=90);
    Modelica_Fluid.Interfaces.FluidPort_ArrayIcon sidePorts[n_sidePorts](redeclare 
        package Medium = Medium, m_flow(each start=0), mXi_flow(each start=0)) 
      annotation (extent=[100,30; 108,-30]);
    
    Medium.BaseProperties medium(
      preferredMediumStates=true,
      p(start=p_ambient),
      T(start=T_start),
      h(start=h_start),
      Xi(start=X_start[1:Medium.nXi]));
    
    SI.Height level(stateSelect=StateSelect.prefer, start=level_start) 
      "Level height of tank";
    SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";
    SI.Energy U "Internal energy of tank volume";
    SI.Mass m "Mass of fluid in tank";
    SI.Mass mXi[Medium.nXi] "Masses of independent components in the fluid";
    
  protected 
    Real Q_lost = 0 "Wärmeverlust (currently dummy)";
    
    outer Modelica_Fluid.Components.FluidOptions fluidOptions 
      "Global default options";
    parameter Medium.SpecificEnthalpy h_start = Medium.specificEnthalpy_pTX(p_ambient, T_start, X_start);
    parameter Types.Init.Temp initOption2=
        if initOption == Types.InitWithGlobalDefault.UseGlobalFluidOption then 
             fluidOptions.default_initOption else initOption 
        annotation(Evaluate=true, Hide=true);
    parameter Integer precision = 3 "Precision for tank level in animation" annotation(Hide=false);
    
    Medium.EnthalpyFlowRate H_flow_topPorts[n_topPorts];
    Medium.EnthalpyFlowRate H_flow_bottomPorts[n_bottomPorts];
    Medium.EnthalpyFlowRate H_flow_sidePorts[n_sidePorts];
    
    Medium.MassFlowRate m_flow_topPorts[n_topPorts];
    Medium.MassFlowRate m_flow_bottomPorts[n_bottomPorts];
    Medium.MassFlowRate m_flow_sidePorts[n_sidePorts];
    
    Medium.MassFlowRate mXi_flow_topPorts[n_topPorts,Medium.nXi];
    Medium.MassFlowRate mXi_flow_bottomPorts[n_bottomPorts,Medium.nXi];
    Medium.MassFlowRate mXi_flow_sidePorts[n_sidePorts,Medium.nXi];
    
  protected 
    Modelica_Fluid.Utilities.TankAttachment tankAttachmentTop[n_topPorts](
      each h=medium.h,
      each d=medium.d,
      each Xi=medium.Xi,
      each p_ambient=p_ambient,
      each level=level,
      each h_start = h_start,
      each X_start = X_start,
      each level_start = level_start,
      each onlyInFlow=true,
      H_flow=H_flow_topPorts,
      m_flow=m_flow_topPorts,
      mXi_flow=mXi_flow_topPorts,
      pipeHeight=top_heights,
      redeclare package Medium = Medium) 
        annotation (extent=[-20,80; 20,40]);
    Modelica_Fluid.Utilities.TankAttachment tankAttachmentBottom[n_bottomPorts](
      each h=medium.h,
      each d=medium.d,
      each Xi=medium.Xi,
      each p_ambient=p_ambient,
      each level=level,
      each h_start = h_start,
      each X_start = X_start,
      each level_start = level_start,
      each onlyInFlow=false,
      H_flow=H_flow_bottomPorts,
      m_flow=m_flow_bottomPorts,
      mXi_flow=mXi_flow_bottomPorts,
      pipeHeight=bottom_heights,
      pipeDiameter=bottom_diameters,
      redeclare package Medium = Medium) 
        annotation (extent=[-20,-80; 20,-40]);
    Modelica_Fluid.Utilities.TankAttachment tankAttachmentSide[n_sidePorts](
      each h=medium.h,
      each d=medium.d,
      each Xi=medium.Xi,
      each p_ambient=p_ambient,
      each level=level,
      each h_start = h_start,
      each X_start = X_start,
      each level_start = level_start,
      each onlyInFlow=false,
      H_flow=H_flow_sidePorts,
      m_flow=m_flow_sidePorts,
      mXi_flow=mXi_flow_sidePorts,
      pipeHeight=side_heights,
      pipeDiameter=side_diameters,
      redeclare package Medium = Medium) 
        annotation (extent=[40,-20; 80,20], rotation=90);
    
  equation 
    for i in 1:n_bottomPorts loop
      connect(tankAttachmentBottom[i].port, bottomPorts[i]) 
        annotation (points=[0,-80.4; 0,-104],
                                            style(color=3, rgbcolor={0,0,255}));
    end for;
    
    for i in 1:n_topPorts loop
       connect(tankAttachmentTop[i].port, topPorts[i]) 
         annotation (points=[0,80.4; 0,104],
                                           style(color=3, rgbcolor={0,0,255}));
    end for;
    
    for i in 1:n_sidePorts loop
       connect(tankAttachmentSide[i].port, sidePorts[i]) 
         annotation (points=[80.4,-1.2491e-015; 104,0],
                                             style(color=3, rgbcolor={0,0,255}));
    end for;
    
    medium.p = p_ambient;
    V = area*level+V0 "Volume of fluid";
    m = V*medium.d "Mass of fluid";
    mXi = m*medium.Xi "Mass of fluid components";
    U = m*medium.u "Internal energy of fluid";
    // Q_lost = - k*2*sqrt(Modelica.Constants.pi*area)*level*(medium.T - T_ambient);
    
    // Mass balances
    der(m) = sum(m_flow_bottomPorts) + sum(m_flow_sidePorts) + sum(m_flow_topPorts);
    for i in 1:Medium.nXi loop
         der(mXi[i]) = sum(mXi_flow_bottomPorts[i,:]) +
                       sum(mXi_flow_sidePorts[i,:]) +
                       sum(mXi_flow_topPorts[i,:]);
    end for;
    
    // Energy balance
    if Medium.singleState then
       der(U) = sum(H_flow_bottomPorts)+sum(H_flow_sidePorts)+sum(H_flow_topPorts) + Q_lost 
        "Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)";
    else
       der(U) = sum(H_flow_bottomPorts)+sum(H_flow_sidePorts)+sum(H_flow_topPorts) - p_ambient*der(V) + Q_lost;
    end if;
    assert(level <= height, "Tank is full (level = height = " + String(level) + ")");
    
  initial equation 
    if initOption2 == Types.Init.NoInit then
      // no initial equations
    elseif initOption2 == Types.Init.InitialValues then
      level = level_start;
      medium.T = T_start;
      medium.Xi = X_start[1:Medium.nXi];
    elseif initOption2 == Types.Init.SteadyState then
      der(level) = 0;
      der(medium.h) = 0;
      der(medium.Xi) = zeros(Medium.nXi);
    elseif initOption2 == Types.Init.SteadyStateHydraulic then
      der(level) = 0;
      medium.T = T_start;
      medium.Xi = X_start[1:Medium.nXi];
    else
      assert(false, "Unsupported initialization option");
    end if;
    annotation (defaultComponentName="tank",
      Icon(
        Rectangle(extent=[-100,100; 100,0],     style(color=7, fillColor=7)),
        Rectangle(extent=DynamicSelect([-100,-100; 100,10],
                                       [-100, -100; 100, (-100 + 200*level/height)]),
            style(
            color=69,
            rgbcolor={0,127,255},
            fillColor=71,
            rgbfillColor={85,170,255},
            fillPattern=1)),
        Line(points=[-100, 100; -100, -100; 100, -100; 100, 100], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=69,
            rgbfillColor={0,127,255},
            fillPattern=1)),
        Text(
          extent=[-95,90; 95,60],
          string="%name",
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=69,
            rgbfillColor={0,127,255},
            fillPattern=1)),
        Text(
          extent=[-95,50; 95,35],
          style(color=0),
          string="%level_start",
        Line(points=[-100,100; 100,100], style(
            color=0,
            rgbcolor={0,0,0},
            pattern=3))),
        Text(
          extent=[-95,30; 95,5],
          style(color=0),
          string=DynamicSelect(" ",realString(level,1,integer(precision)))),
        Line(points=[-100,100; 100,100], style(
            color=0,
            rgbcolor={0,0,0},
            pattern=3)),
        Text(
          extent=[-40,124; -22,112],
          string="1",
          style(color=10, rgbcolor={135,135,135})),
        Text(
          extent=[111,37; 129,25],
          string="1",
          style(color=10, rgbcolor={135,135,135})),
        Text(
          extent=[-38,-112; -20,-124],
          string="1",
          style(color=10, rgbcolor={135,135,135}))),
      Documentation(info="<HTML>
<p>
This is a simplified model of a tank. 
The top part is open to the environment at the fixed pressure 
<tt>p_ambient</tt>. Heat transfer to the environment and to 
the tank walls is neglected.
The tank is filled with a single or multiple-substance liquid, 
assumed to have uniform temperature and mass fractions.
</p>
 
<p>
There are three arrays of FluidConnectors (topPorts, bottomPorts, sidePorts)
at which other fluid components can be attached. The difference between
these ports is just the graphical layout on the icon, to visualize
whether a pipe is connected to the top, bottom or side of the tank.
All ports are defined by height (with respect to bottom of tank)
and the diameter of the component that is attached to this port.
</p>
 
<p>
In a diagram animation, the fill level of the tank
is visualized, as well as the value of level.
</p>
 
</HTML>", revisions="<html>
<ul>
<li><i>Jan. 6, 2006</i>
   implemented by Katja Poschlad, Manuel Remelhel (AST Uni Dortmund), 
   Martin Otter (DLR), based on Components.Tank model.</li>
</ul>
</html>"),
      Diagram,
      uses(Modelica(version="2.2.1"), Modelica_Fluid(version="0.952")),
      Coordsys(grid=[1,1], scale=0.2));
  end OpenTank;
  
  model ValveIncompressible "Valve for (almost) incompressible fluids" 
    extends Interfaces.PartialValve;
    import Modelica_Fluid.Types.CvTypes;
  annotation (
    Icon(Text(extent=[-100, -40; 100, -80], string="%name")),
    Diagram,
    Documentation(info="<HTML>
<p>Valve model according to the IEC 534/ISA S.75 standards for valve sizing, incompressible fluids. <p>
Extends the <tt>Interfaces.PartialValve</tt> model (see the corresponding documentation for common valve features).
<p>This model can be used with any low compressibility fluids, such as liquids or gases at very low pressure drops.
</html>",
      revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"),
      Coordsys(grid=[1,1], scale=0));
  initial equation 
    if CvData == CvTypes.OpPoint then
      m_flow_nom = flowCharacteristic(stemPosition_nom)*Av*sqrt(d_nom)*sqrtR(dp_nom) 
        "Determination of Av by the operating point";
    end if;
    
  equation 
    if CheckValve then
      m_flow = flowCharacteristic(stemPosition)*Av*sqrt(d)*
               smooth(0,if dp>=0 then sqrtR(dp) else 0);
    else
      m_flow = flowCharacteristic(stemPosition)*Av*sqrt(d)*sqrtR(dp);
    end if;
  end ValveIncompressible;
  
  model ValveVaporizing 
    "Valve for possibly vaporizing (almost) incompressible fluids, accounts for choked flow conditions" 
    import Modelica_Fluid.Types.CvTypes;
    extends Interfaces.PartialValve(
      redeclare replaceable package Medium = 
      Modelica.Media.Interfaces.PartialTwoPhaseMedium);
    parameter Real Fl_nom=0.9 "Liquid pressure recovery factor";
    replaceable function FlCharacteristic = 
        Modelica_Fluid.Types.ValveCharacteristics.one 
      extends Modelica_Fluid.Types.ValveCharacteristics.baseFun 
      "Pressure recovery characteristic";
    Real Ff "Ff coefficient (see IEC/ISA standard)";
    Real Fl "Pressure recovery coefficient Fl (see IEC/ISA standard)";
    Medium.AbsolutePressure pv "Saturation pressure";
    SI.Pressure dpEff "Effective pressure drop";
    Medium.AbsolutePressure pin "Inlet pressure";
    Medium.AbsolutePressure pout "Outlet pressure";
    annotation (
      Icon(Text(extent=[-100, -40; 100, -80], string="%name")),
      Diagram,
      Documentation(info="<HTML>
<p>Valve model according to the IEC 534/ISA S.75 standards for valve sizing, incompressible fluid, with possible choked flow conditions. <p>
Extends the <tt>Interfaces.PartialValve</tt> model (see the corresponding documentation for common valve features).<p>
The model operating range includes choked flow operation, which takes place for low outlet pressures due to flashing in the vena contracta; otherwise, non-choking conditions are assumed.
<p>This model can be used with two-phase medium models, to describe the liquid and (possible) two-phase conditions.
<p>The default liquid pressure recovery coefficient <tt>Fl</tt> is constant and given by the parameter <tt>Fl_nom</tt>. The relative change (per unit) of the recovery coefficient can be specified as a given function of the valve opening by replacing the <tt>FlCharacteristic</tt> function.
</HTML>",
        revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"));
  initial equation 
    if CvData == CvTypes.OpPoint then
      m_flow_nom = flowCharacteristic(stemPosition_nom)*Av*sqrt(d_nom)*sqrtR(dp_nom) 
        "Determination of Av by the operating point";
    end if;
  equation 
    pin = port_a.p;
    pout = port_b.p;
    pv = Medium.saturationPressure(T);
    Ff = 0.96 - 0.28*sqrt(pv/Medium.fluidConstants[1].criticalPressure);
    Fl = Fl_nom*FlCharacteristic(stemPosition);
    dpEff = if pout < (1 - Fl^2)*pin + Ff*Fl^2*pv then 
              Fl^2*(pin - Ff*pv) else dp 
      "Effective pressure drop, accounting for possible choked conditions";
    if CheckValve then
       m_flow = flowCharacteristic(stemPosition)*Av*sqrt(d)*
           (if dpEff>=0 then sqrtR(dpEff) else 0);
     else
       m_flow = flowCharacteristic(stemPosition)*Av*sqrt(d)*sqrtR(dpEff);
    end if;
  end ValveVaporizing;
  
  model ValveCompressible 
    "Valve for compressible fluids, accounts for choked flow conditions" 
    extends Interfaces.PartialValve;
    import Modelica_Fluid.Types.CvTypes;
    parameter Real Fxt_full=0.5 "Fk*xt critical ratio at full opening";
    replaceable function xtCharacteristic = 
        Modelica_Fluid.Types.ValveCharacteristics.one 
      extends Modelica_Fluid.Types.ValveCharacteristics.baseFun 
      "Critical ratio characteristic";
    Real Fxt;
    Real x "Pressure drop ratio";
    Real xs "Saturated pressure drop ratio";
    Real Y "Compressibility factor";
    Medium.AbsolutePressure p "Inlet pressure";
  protected 
    parameter Real Fxt_nom(fixed=false) "Nominal Fxt";
    parameter Real x_nom(fixed=false) "Nominal pressure drop ratio";
    parameter Real xs_nom(fixed=false) "Nominal saturated pressure drop ratio";
    parameter Real Y_nom(fixed=false) "Nominal compressibility factor";
    
    annotation (
    Icon(Text(extent=[-100, -40; 100, -80], string="%name")),
    Diagram,
    Documentation(info="<HTML>
<p>Valve model according to the IEC 534/ISA S.75 standards for valve sizing, compressible fluid. <p>
Extends the <tt>Interfaces.PartialValve</tt> model (see the corresponding documentation for common valve features).
<p>The product Fk*xt is given by the parameter <tt>Fxt_full</tt>, and is assumed constant by default. The relative change (per unit) of the xt coefficient with the valve opening can be specified by replacing the <tt>xtCharacteristic</tt> function.
</HTML>",
      revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"));
  initial equation 
    if CvData == CvTypes.OpPoint then
      // Determination of Av by the nominal operating point conditions
      Fxt_nom = Fxt_full*xtCharacteristic(stemPosition_nom);
      x_nom = dp_nom/p_nom;
      xs_nom = smooth(0, if x_nom > Fxt_nom then Fxt_nom else x_nom);
      Y_nom = 1 - abs(xs_nom)/(3*Fxt_nom);
      m_flow_nom = flowCharacteristic(stemPosition_nom)*Av*Y_nom*sqrt(d_nom)*sqrtR(p_nom*xs_nom);
    else
      // Dummy values
      Fxt_nom = 0;
      x_nom = 0;
      xs_nom = 0;
      Y_nom = 0;
    end if;
    
  equation 
    p = noEvent(if dp>=0 then port_a.p else port_b.p);
    Fxt = Fxt_full*xtCharacteristic(stemPosition);
    x = dp/p;
    xs = smooth(0, if x < -Fxt then -Fxt else if x > Fxt then Fxt else x);
    Y = 1 - abs(xs)/(3*Fxt);
    if CheckValve then
      m_flow = flowCharacteristic(stemPosition)*Av*Y*sqrt(d)*
        smooth(0,if xs>=0 then sqrtR(p*xs) else 0);
    else
      m_flow = flowCharacteristic(stemPosition)*Av*Y*sqrt(d)*sqrtR(p*xs);
    end if;
  end ValveCompressible;
  
  model ValveLinear "Valve for water/steam flows with linear pressure drop" 
    extends Interfaces.PartialTwoPortTransport;
    parameter Types.HydraulicConductance Kv 
      "Hydraulic conductance at full opening";
    Modelica.Blocks.Interfaces.RealInput opening 
    annotation (extent=[-20,70; 20,110],   rotation=-90);
  equation 
    m_flow = Kv*opening*dp;
    
  annotation (
    Icon(
        Polygon(points=[-100,50; -100,-50; 0,0; -100,50],  style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Line(points=[0,60; 0,0],   style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Rectangle(extent=[-20,70; 20,50],   style(
            color=0,
            fillColor=0,
            fillPattern=1)),
        Polygon(points=[100,50; 0,0; 100,-50; 100,50],  style(
            color=0,
            thickness=2,
            fillPattern=1)),
           Text(extent=[-143,-66; 148,-106],  string="%name")),
    Diagram,
    Documentation(info="<HTML>
<p>This very simple model provides a pressure drop which is proportional to the flowrate and to the <tt>opening</tt> signal, without computing any fluid property.
<p>A medium model must be nevertheless be specified, so that the fluid ports can be connected to other components using the same medium model.
</HTML>",
      revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"),
      Coordsys(grid=[1,1], scale=0));
  end ValveLinear;
  
  model ValveDiscrete "Valve for water/steam flows with linear pressure drop" 
    extends Modelica_Fluid.Interfaces.PartialTwoPortTransport;
    parameter Modelica_Fluid.Types.HydraulicConductance Kv 
      "Hydraulic conductance for open valve (m_flow = Kv*dp)";
    parameter Real Kv_small_rel = 0 
      "Relative hydraulic conductance for closed valve (m_flow = Kv_small_rel*Kv*dp)";
    Modelica.Blocks.Interfaces.BooleanInput open 
    annotation (extent=[-20,60; 20,100],   rotation=-90);
  equation 
    m_flow = if open then Kv*dp else Kv_small_rel*Kv*dp;
    
  annotation (
    Icon(
        Line(points=[0,50; 0,0],   style(
            color=0,
            rgbcolor={0,0,0},
            fillPattern=1)),
        Rectangle(extent=[-20,60; 20,50],   style(
            color=0,
            fillColor=0,
            fillPattern=1)),
           Text(extent=[-145,-58; 146,-98],   string="%name"),
        Polygon(points=[-100,50; 100,-50; 100,50; 0,0; -100,-50; -100,50], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=DynamicSelect(7, if open > 0.5 then 2 else 7)))),
    Diagram,
    Documentation(info="<HTML>
<
<p>
This very simple model provides a pressure drop which is proportional to the flowrate if the Boolean open signal is <b>true</b>. Otherwise, the
mass flow rate is zero. If Kv_small_rel > 0, a small leakage
mass flow rate occurs when open = <b>false</b>. This might be
useful in certain situations when the model is not
mathematically well-defined due to a closed valve.
</p>

<p>
In a diagram animation, the valve is shown in \"green\", when
it is open.
</p>

</HTML>",
      revisions="<html>
<ul>
<li><i>Nov 2005</i>
    by Katja Poschlad (based on ValveLinear).</li>
</ul>
</html>"),
      Coordsys(grid=[1,1], scale=0));
  end ValveDiscrete;
  
  model Pump "Centrifugal pump with ideally controlled speed" 
    extends Interfaces.PartialPump;
    import Modelica.SIunits.Conversions.NonSIunits.*;
    parameter AngularVelocity_rpm N_const = N_nom "Constant rotational speed";
    Modelica.Blocks.Interfaces.RealInput N_in(redeclare type SignalType = 
          Modelica.SIunits.Conversions.NonSIunits.AngularVelocity_rpm) 
      "Prescribed rotational speed" 
      annotation (extent=[-36,34; -16,54],   rotation=-90);
  equation 
    N = N_in "Rotational speed";
    if cardinality(N_in)==0 then
      N_in = N_const "Rotational speed provided by parameter";
    end if;
    annotation (
      Icon(
        Text(extent=[-58,58; -30,38], string="n")),
      Diagram,
      Documentation(info="<HTML>
<p>This model describes a centrifugal pump (or a group of <tt>Np</tt> pumps in parallel) with controlled speed, either fixed or provided by an external signal.
<p>The model extends <tt>PartialPump</tt>
<p>If the <tt>N_in</tt> input connector is wired, it provides rotational speed of the pumps (rpm); otherwise, a constant rotational speed equal to <tt>n_const</tt> (which can be different from <tt>N_nom</tt>) is assumed.</p>
</HTML>",
        revisions="<html>
<ul>
<li><i>31 Oct 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  end Pump;
  
  model PumpShaft "Centrifugal pump with mechanical connector for the shaft" 
    extends Interfaces.PartialPump;
    SI.Angle phi "Shaft angle";
    SI.AngularVelocity omega "Shaft angular velocity";
    Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft 
    annotation (extent=[80,4; 110,32]);
  equation 
    phi = shaft.phi;
    omega = der(phi);
    N = Modelica.SIunits.Conversions.to_rpm(omega);
    W_single = omega*shaft.tau;
  annotation (
    Icon(Rectangle(extent=[60,26; 84,12], style(
            color=10,
            rgbcolor={95,95,95},
            gradient=2,
            fillColor=10,
            rgbfillColor={95,95,95}))),
    Diagram,
    Documentation(info="<HTML>
<p>This model describes a centrifugal pump (or a group of <tt>Np</tt> pumps in parallel) with a mechanical rotational connector for the shaft, to be used when the pump drive has to be modelled explicitly. In the case of <tt>Np</tt> pumps in parallel, the mechanical connector is relative to a single pump.
<p>The model extends <tt>PartialPump</tt>
 </HTML>",
       revisions="<html>
<ul>
<li><i>31 Oct 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  end PumpShaft;
  
  model Evaporator 
    "Simple Evaporator with two states, see Astroem, Bell: Drum-boiler dynamics, Automatica 36, 2000, pp.363-378" 
    import Modelica.SIunits.Conversions.*;
    import Modelica.Constants;
    import Modelica_Fluid.Types;
    import Modelica_Fluid.Types.FlowDirection;
    import Modelica_Fluid.Types.FlowDirectionWithGlobalDefault;
    replaceable package Medium = 
        Modelica.Media.Interfaces.PartialTwoPhaseMedium 
      extends Modelica.Media.Interfaces.PartialTwoPhaseMedium "Medium model" 
        annotation (choicesAllMatching=true);
    parameter SI.Mass m_D "mass of surrounding drum metal";
    parameter Medium.SpecificHeatCapacity cp_D 
      "specific heat capacity of drum metal";
    parameter SI.Volume V_t "total volume inside drum";
    parameter Types.InitWithGlobalDefault.Temp initOption=
              Types.InitWithGlobalDefault.UseGlobalFluidOption 
      "Initialization option" 
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_start = Medium.p_default 
      "Start value of pressure" 
      annotation(Dialog(tab = "Initialization"));
    parameter SI.Volume V_l_start = V_t/2 
      "Start value of liquid volumeStart value of volume" 
      annotation(Dialog(tab = "Initialization"));
    
    parameter FlowDirectionWithGlobalDefault.Temp flowDirection=
              FlowDirectionWithGlobalDefault.UseGlobalFluidOption 
      "Unidirectional (port_a -> port_b) or bidirectional flow component" 
       annotation(Dialog(tab="Advanced"));
    
    Interfaces.FluidPort_a feedwater(redeclare package Medium = Medium,
                       m_flow(min=if allowFlowReversal then -Constants.inf else 0)) 
      annotation (extent=[-110,-10; -90,10]);
    Interfaces.FluidPort_b steam(redeclare package Medium = Medium,
                       m_flow(max=if allowFlowReversal then +Constants.inf else 0)) 
      annotation (extent=[110,-10; 90,10]);
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort 
      annotation (extent=[-10,-110; 10,-90]);
    Modelica.Blocks.Interfaces.RealOutput V(
      redeclare type SignalType = SI.Volume) "liquid volume" 
      annotation (extent=[30, 100; 50, 120], rotation=90);
    
    Medium.SaturationProperties sat 
      "State vector to compute saturation properties";
    Medium.AbsolutePressure p(start=p_start, stateSelect=StateSelect.prefer) 
      "pressure inside drum boiler";
    Medium.Temperature T "temperature inside drum boiler";
    SI.Volume V_v "volume of vapour phase";
    SI.Volume V_l(start=V_l_start, stateSelect=StateSelect.prefer) 
      "volumes of liquid phase";
    Medium.SpecificEnthalpy h_v=Medium.dewEnthalpy(sat) 
      "specific enthalpy of vapour";
    Medium.SpecificEnthalpy h_l=Medium.bubbleEnthalpy(sat) 
      "specific enthalpy of liquid";
    Medium.Density rho_v=Medium.dewDensity(sat) "density in vapour phase";
    Medium.Density rho_l=Medium.bubbleDensity(sat) "density in liquid phase";
    SI.Mass m "total mass of drum boiler";
    SI.Energy U "internal energy";
    Medium.Temperature T_D=heatPort.T "temperature of drum";
    SI.HeatFlowRate q_F=heatPort.Q_flow "heat flow rate from furnace";
    Medium.SpecificEnthalpy h_W=feedwater.h "feed water enthalpy";
    Medium.SpecificEnthalpy h_S=steam.h "steam enthalpy";
    SI.MassFlowRate qm_W=feedwater.m_flow "feed water mass flow rate";
    SI.MassFlowRate qm_S=steam.m_flow "steam mass flow rate";
  protected 
    outer Modelica_Fluid.Components.FluidOptions fluidOptions 
      "Global default options";
    parameter Types.Init.Temp initOption2=
        if initOption == Types.InitWithGlobalDefault.UseGlobalFluidOption then 
             fluidOptions.default_initOption else initOption 
        annotation(Evaluate=true, Hide=true);
    parameter Boolean allowFlowReversal=
       flowDirection == FlowDirectionWithGlobalDefault.Bidirectional
       or flowDirection == FlowDirectionWithGlobalDefault.UseGlobalFluidOption
       and fluidOptions.default_flowDirection ==FlowDirection.Bidirectional 
      "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
       annotation(Evaluate=true, Hide=true);
  equation 
    // balance equations  
    m = rho_v*V_v + rho_l*V_l + m_D "Total mass";
    U = rho_v*V_v*h_v + rho_l*V_l*h_l - p*V_t + m_D*cp_D*T_D "Total energy";
    der(m) = qm_W + qm_S "Mass balance";
    der(U) = q_F + qm_W*h_W + qm_S*h_S "Energy balance";
    V_t = V_l + V_v;
    
    // Properties of saturated liquid and steam
    sat.psat = p;
    sat.Tsat = T;
    sat.Tsat = Medium.saturationTemperature(p);
    
    // ideal heat transfer between metal and water
    T_D = T;
    
    // boundary conditions at the ports
    feedwater.p = p;
    feedwater.H_flow = semiLinear(feedwater.m_flow, feedwater.h, h_l);
    steam.p = p;
    steam.H_flow = semiLinear(steam.m_flow, steam.h, h_v);
    
    // liquid volume 
    V = V_l;
    
    // Check that two-phase equilibrium is actually possible
    assert(p<Medium.fluidConstants[1].criticalPressure-10000,
           "Evaporator model requires subcritical pressure");
  initial equation 
    // Initial conditions
    if initOption2 == Types.Init.NoInit then
      // no initial equations
    elseif initOption2 == Types.Init.InitialValues then
     p = p_start;
     V_l = V_l_start;
    elseif initOption2 == Types.Init.SteadyState then
      der(p) = 0;
      der(V_l) = 0;
    elseif initOption2 == Types.Init.SteadyStateHydraulic then
      der(p) = 0;
      V_l = V_l_start;
    else
      assert(false, "Unsupported initialization option");
    end if;
    
    annotation (
      Coordsys(grid=[1, 1], component=[20, 20]),
      Diagram,
      Icon(
        Rectangle(extent=[-100, 59; 100, -61], style(
            color=0,
            gradient=2,
            fillColor=8)),
        Rectangle(extent=[-100, 34; 100, -36], style(
            color=69,
            gradient=2,
            fillColor=69)),
        Ellipse(extent=[18, 0; 48, -29], style(pattern=0, fillColor=7)),
        Ellipse(extent=[-1, 29; 29, 0], style(pattern=0, fillColor=7)),
        Ellipse(extent=[48, 34; 78, 5], style(pattern=0, fillColor=7)),
        Ellipse(extent=[-31, 1; -1, -28], style(pattern=0, fillColor=7)),
        Ellipse(extent=[47, 14; 77, -15], style(pattern=0, fillColor=7)),
        Ellipse(extent=[-72, 25; -42, -4], style(pattern=0, fillColor=7)),
        Ellipse(extent=[71, 0; 101, -29], style(pattern=0, fillColor=7)),
        Ellipse(extent=[74, 14; 104, -15], style(pattern=0, fillColor=7)),
        Ellipse(extent=[71, 29; 101, 0], style(pattern=0, fillColor=7)),
        Text(
          extent=[-139,111; 144,57],
          string="%name",
          style(gradient=2, fillColor=69)),
        Line(points=[0, -60; 0, -100], style(color=42)),
        Line(points=[40, 99; 40, 60])),
      Documentation(revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
     Initialization options fixed</li>
<li><i>6 Sep 2005</i><br>
    Model by Ruediger Franke modified after the 45th Design Meeting</li>
</ul>
</html>", info="<html>
Model of a simple evaporator with two states. The model assumes two-phase equilibrium inside the component; saturated steam goes out of the steam outlet.
<p>
References: Astroem, Bell: Drum-boiler dynamics, Automatica 36, 2000, pp.363-378
</html>"));
  equation 
    
  end Evaporator;
  
  model Tank "Obsolet component (use instead Components.OpenTank)" 
    import Modelica_Fluid.Types;
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
      annotation (choicesAllMatching=true);
    parameter SI.Area area "Tank area";
    parameter SI.Area pipeArea "Area of outlet pipe";
    parameter SI.Volume V0 = 0 "Volume of the liquid when the level is zero";
    parameter SI.Height H0 = 0 
      "Height of zero level reference over the bottom port";
    parameter Medium.AbsolutePressure p_ambient=fluidOptions.default_p_ambient 
      "Tank surface pressure";
    parameter Types.InitWithGlobalDefault.Temp initOption=
              Types.InitWithGlobalDefault.UseGlobalFluidOption 
      "Initialization option" 
      annotation(Dialog(tab = "Initialization"));
    parameter Boolean use_T_start = true 
      "Use T_start if true, otherwise h_start" 
      annotation(Dialog(tab = "Initialization"), Evaluate = true);
    parameter Medium.Temperature T_start=
      if use_T_start then Medium.T_default else Medium.temperature_phX(p_ambient,h_start,X_start) 
      "Start value of temperature" 
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.specificEnthalpy_pTX(p_ambient, T_start, X_start) else Medium.h_default 
      "Start value of specific enthalpy" 
      annotation(Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    parameter SI.Height level_start(min=0) "Initial tank level" 
      annotation(Dialog(tab="Initialization"));
    
    Interfaces.FluidPort_b port(redeclare package Medium = Medium,
                                m_flow(start=0), mXi_flow(each start=0)) 
      annotation (extent=[-10, -120; 10, -100], rotation=90);
    Medium.BaseProperties medium(
      preferredMediumStates=true,
      p(start=p_ambient),
      T(start=T_start),
      Xi(start=X_start[1:Medium.nXi]));
    
    SI.Height level(start=level_start,stateSelect=StateSelect.prefer) 
      "Height of tank level over the zero reference";
    Medium.AbsolutePressure p_bottom "Pressure at bottom of tank";
    SI.Energy U "Internal energy of tank volume";
    SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";
    Real m(quantity=Medium.mediumName, unit="kg", stateSelect=StateSelect.never) 
      "Mass of tank volume";
    Real mX[Medium.nX](quantity=Medium.substanceNames, each unit="kg") 
      "Component masses of the independent substances";
  protected 
    outer Modelica_Fluid.Components.FluidOptions fluidOptions 
      "Global default options";
    parameter Types.Init.Temp initOption2=
        if initOption == Types.InitWithGlobalDefault.UseGlobalFluidOption then 
             fluidOptions.default_initOption else initOption 
        annotation(Evaluate=true, Hide=true);
  equation 
    medium.p = p_ambient;
    V = area*level+V0 "Volume of fluid";
    m = V*medium.d "Mass of fluid";
    mX = m*medium.Xi "Mass of fluid components";
    U = m*medium.u "Internal energy of fluid";
    
    // Mass balance
    der(m) = port.m_flow;
    der(mX) = port.mXi_flow;
    
    // Momentum balance
    p_bottom = (medium.d*fluidOptions.g*(level+H0)) + p_ambient;
    port.p = p_bottom - smooth(2,noEvent(if port.m_flow < 0 then port.m_flow^2/(2*medium.d*pipeArea^2) else 0));
    
    // Energy balance
    if Medium.singleState then
      der(U) = port.H_flow "Mechanical work is neglected";
    else
      der(U) = port.H_flow - p_ambient*der(V);
    end if;
    
    /* Handle reverse and zero flow */
    port.H_flow = semiLinear(port.m_flow, port.h, medium.h);
    port.mXi_flow = semiLinear(port.m_flow, port.Xi, medium.X);
    
  initial equation 
    if initOption2 == Types.Init.NoInit then
      // no initial equations
    elseif initOption2 == Types.Init.InitialValues then
      level = level_start;
      if use_T_start then
        medium.T = T_start;
      else
        medium.h = h_start;
      end if;
      medium.Xi = X_start[1:Medium.nXi];
    elseif initOption2 == Types.Init.SteadyState then
      der(level) = 0;
      der(medium.h) = 0;
      der(medium.Xi) = zeros(Medium.nXi);
    elseif initOption2 == Types.Init.SteadyStateHydraulic then
      der(level) = 0;
      if use_T_start then
        medium.T = T_start;
      else
        medium.h = h_start;
      end if;
      medium.Xi = X_start[1:Medium.nXi];
    else
      assert(false,"Unsupported initialization option initOption = " + String(initOption2)
                   +"\nin model Modelica_Fluid.Components.Tank");
    end if;
    annotation (
      Icon(
        Rectangle(extent=[-100, 90; 100, 26], style(color=7, fillColor=7)),
        Rectangle(extent=[-100, 26; 100, -100], style(
            color=69,
            fillColor=69,
            fillPattern=1)),
        Line(points=[-100, 100; -100, -100; 100, -100; 100, 100], style(
            color=0,
            fillColor=69,
            fillPattern=1)),
        Text(
          extent=[-112, 162; 122, 102],
          string="%name",
          style(fillColor=69, fillPattern=1)),
        Text(
          extent=[-86, -38; 94, -78],
          style(color=0),
          string="%level_start"),
        Text(
          extent=[-94, 78; 94, 38],
          style(color=0),
          string="%p_ambient"),
        Text(
          extent=[-94, 14; 90, -2],
          style(color=0),
          string="level_start")),
      Documentation(info="<HTML>
<p>
This is a simplified model of a tank. The top part is open to the environment at the fixed pressure <tt>p_ambient</tt>. Heat transfer to the environment and to the tank walls is neglected.
The tank is filled with a single or multiple-substance liquid, assumed to have uniform temperature and mass fractions.
<p>The geometry of the tank is specified by the following parameters: <tt>V0</tt> is the volume of the liquid when the level is at the zero reference; <tt>area</tt> is the cross-sectional area of the tank; <tt>H0</tt> is the height of the zero-reference level plane over the port connector. It is thus possible to model rounded-bottom tanks, as long as they have a cylindrical shape in the range of operating levels.
<p>The tank can be initialized with the following options:
<ul>
<li>NoInit: no explicit initial conditions
<li>InitialValues: initial values of temperature (or specific enthalpy), composition and level are specified
<li>SteadyStateHydraulic: initial values of temperature (or specific enthalpy) and composition are specified; the initial level is determined so that levels and pressure are at steady state.
</ul>
Full steady state initialization is not supported, because the corresponding intial equations for temperature/enthalpy are undetermined (the flow rate through the port at steady state is zero). 
</p>
</HTML>",   revisions="<html>
<ul>
<li><i>1 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from a previous version of Modelica_Fluid</li>
</ul>
</html>"),
      Diagram);
  end Tank;
  
  model PressureDropPipe 
    "Obsolet component (use instead PressureLosses.WallFrictionAndGravity)" 
    extends Modelica_Fluid.Interfaces.PartialTwoPortTransport;
    extends Modelica_Fluid.Utilities.PipeFriction;
    annotation (Icon(
        Rectangle(extent=[-100,60; 100,-60],   style(
            color=0,
            gradient=2,
            fillColor=8)),
        Rectangle(extent=[-100,34; 100,-36],   style(
            color=69,
            gradient=2,
            fillColor=69)),
        Text(
          extent=[-150,140; 150,80],
          string="%name",
          style(gradient=2, fillColor=69))), Documentation(info="<html>
<p>
This model describes pressure losses due to friction in a pipe. It is assumed that no mass or energy is stored in the pipe. 
The details of the pipe friction model are described
<a href=\"Modelica://Modelica_Fluid.Utilities.PipeFriction\">here</a>.
</p>
</html>"));
  equation 
    if frictionType == Modelica_Fluid.Types.FrictionTypes.DetailedFriction then
       if from_dp then
          d = if dp >= 0 then medium_a.d else medium_b.d;
          eta = if dp >= 0 then Medium.dynamicViscosity(medium_a) else 
                               Medium.dynamicViscosity(medium_b);
       else
          d = if m_flow >= 0 then medium_a.d else medium_b.d;
          eta = if m_flow >= 0 then Medium.dynamicViscosity(medium_a) else 
                               Medium.dynamicViscosity(medium_b);
       end if;
    else
      // Assign dummy values for auxiliary variables
       d = 0;
       eta = 0;
    end if;
  end PressureDropPipe;
  
model LongPipe "Distributed pipe model with optional wall" 
    import Modelica_Fluid;
    
  extends Modelica_Fluid.Interfaces.Flow1D(
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
    
  replaceable 
      Modelica_Fluid.WorkInProgress.Utilities.PipeHeatTransfer.PipeHT_constAlpha
      heat(
    redeclare final package Medium = Medium,
    final n=n,
    final d_h=d_h,
    final A_h=area_h,
    T=medium.T) extends Modelica_Fluid.HeatTransfer.PartialPipeHeatTransfer(
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
        Modelica_Fluid.Components.Wall_constProps                    extends 
      Modelica_Fluid.Interfaces.PartialPipeWall "Wall model"       annotation(choicesAllMatching, Dialog(enable=use_wall, tab="General", group="Wall - optional"));
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
end LongPipe;
  
model Wall_constProps 
    "Pipe wall, assuming ideal 1D-conduction and constant material properties" 
  extends Modelica_Fluid.Interfaces.PartialPipeWall;
  parameter SI.Density d_wall "Density of wall material";
  parameter SI.SpecificHeatCapacity c_wall 
      "Specific heat capacity of wall material";
  parameter SI.Temperature T_start "Start value for wall temperature";
  parameter SI.Mass[n] m=ones(n)*(a_outer-a_inner)*length*d_wall/n "Wall mass";
  parameter Modelica_Fluid.Types.Init.Temp initOption;
  SI.Temperature[n] T(start=ones(n)*T_start, stateSelect=StateSelect.prefer) 
      "Wall temperature";
initial equation 
  if initOption==3 then
    der(T)=zeros(n);
  else
    T=ones(n)*T_start;
  end if;
equation 
    
  for i in 1:n loop
   assert(m[i]>0, "Wall has negative dimensions");
   c_wall*m[i]*der(T[i]) = thermalPort_a[i].Q_flow + thermalPort_b[i].Q_flow;
  end for;
  //assuming ideal heat conduction perpendicular to fluid flow, conduction in remaining two dimensions is neglected
  thermalPort_a.T=T;
  thermalPort_b.T=T;
end Wall_constProps;
  
model HeatExchanger "Double pipe heat exchanger with outer wall neglected" 
    
  //General
  parameter Integer n(min=1) "Spatial segmentation";
  replaceable package Medium_1 = Modelica.Media.Water.StandardWater extends 
      Modelica.Media.Interfaces.PartialMedium "Inner pipe medium" 
                                                    annotation(choicesAllMatching);
  replaceable package Medium_2 = Modelica.Media.Water.StandardWater extends 
      Modelica.Media.Interfaces.PartialMedium "Outer pipe medium" 
                                                    annotation(choicesAllMatching);
  parameter SI.Length di_1(min=0) "Inner diameter of inner pipe"     annotation(Dialog(tab="General", group="Dimensions"));
  parameter SI.Length da_1(min=0) "Outer diameter of inner pipe"     annotation(Dialog(tab="General", group="Dimensions"));
  parameter SI.Length di_2(min=0) "Inner diameter of outer pipe"     annotation(Dialog(tab="General", group="Dimensions"));
  parameter SI.Length length(min=0) "Length of both pipes" annotation(Dialog(tab="General", group="Dimensions"));
    
  //Wall
  parameter SI.Density d_wall "Density of wall material" annotation(Dialog(tab="General", group="Constant material properties"));
  parameter SI.SpecificHeatCapacity c_wall 
      "Specific heat capacity of wall material"                                      annotation(Dialog(tab="General", group="Constant material properties"));
  final parameter SI.Mass m_wall=sum(pipe_1.wall.m) "Wall mass";
  parameter Boolean initWall_steadyState=false 
      "= true, Wall initialization in steady state"                                    annotation(Dialog(tab="Initialization", group="Wall"));
  parameter SI.Temperature T_start_wall "Start value of wall temperature" annotation(Dialog(tab="Initialization", group="Wall"));
    
  //Initialization pipe 1
  parameter Modelica_Fluid.Types.Init.Temp initOption_1 "Initialization option"
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Inner pipe"));
  parameter Boolean use_T_start_1=true "Use T_start if true, otherwise h_start"
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Inner pipe"));
  parameter Medium_1.AbsolutePressure p_a_start1=Medium_1.p_default 
      "Start value of pressure" 
    annotation(Dialog(tab = "Initialization", group = "Inner pipe"));
  parameter Medium_1.AbsolutePressure p_b_start1=Medium_1.p_default 
      "Start value of pressure" 
    annotation(Dialog(tab = "Initialization", group = "Inner pipe"));
  parameter Medium_1.Temperature T_start_1=if use_T_start_1 then Medium_1.T_default else 
      Medium_1.temperature_phX((p_a_start1+p_b_start1)/2, h_start_1, X_start_1) 
      "Start value of temperature" 
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Inner pipe", enable = use_T_start_1));
  parameter Medium_1.SpecificEnthalpy h_start_1=if use_T_start_1 then 
      Medium_1.specificEnthalpy_pTX((p_a_start1+p_b_start1)/2, T_start_1, X_start_1[1:Medium_1.nXi]) else Medium_1.h_default 
      "Start value of specific enthalpy" 
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Inner pipe", enable = not use_T_start_1));
  parameter Medium_1.MassFraction X_start_1[Medium_1.nX]=Medium_1.X_default 
      "Start value of mass fractions m_i/m" 
    annotation (Dialog(tab="Initialization", group = "Inner pipe", enable=(Medium_1.nXi > 0)));
  parameter Medium_1.MassFlowRate mflow_start_1 "Start value of mass flow rate"
                                    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Inner pipe"));
  //Initialization pipe 2
  parameter Modelica_Fluid.Types.Init.Temp initOption_2 "Initialization option"
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Outer pipe"));
  parameter Boolean use_T_start_2=true "Use T_start if true, otherwise h_start"
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Outer pipe"));
  parameter Medium_2.AbsolutePressure p_a_start2=Medium_2.p_default 
      "Start value of pressure" 
    annotation(Dialog(tab = "Initialization", group = "Outer pipe"));
  parameter Medium_2.AbsolutePressure p_b_start2=Medium_2.p_default 
      "Start value of pressure" 
    annotation(Dialog(tab = "Initialization", group = "Outer pipe"));
  parameter Medium_2.Temperature T_start_2=if use_T_start_2 then Medium_2.T_default else 
      Medium_2.temperature_phX((p_a_start2+p_b_start2)/2, h_start_2, X_start_2) 
      "Start value of temperature" 
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Outer pipe", enable = use_T_start_2));
  parameter Medium_2.SpecificEnthalpy h_start_2=if use_T_start_2 then 
      Medium_2.specificEnthalpy_pTX((p_a_start2+p_b_start2)/2, T_start_2, X_start_2[1:Medium_2.nXi]) else Medium_2.h_default 
      "Start value of specific enthalpy" 
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Outer pipe", enable = not use_T_start_2));
  parameter Medium_2.MassFraction X_start_2[Medium_2.nX]=Medium_2.X_default 
      "Start value of mass fractions m_i/m" 
    annotation (Dialog(tab="Initialization", group = "Outer pipe", enable=Medium_2.nXi>0));
  parameter Medium_2.MassFlowRate mflow_start_2 "Start value of mass flow rate"
                                       annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Outer pipe"));
  //Advanced
  parameter Boolean lumped_dp = true 
      " = true, lumped pressure drop, reduces number of pressure states to one"
                              annotation(Evaluate=true, Dialog(tab="Advanced", group="Conservation of mass, energy, momentum"));
  parameter Boolean static "= true, use quasistatic mass and energy balances" 
                           annotation(Evaluate=true, Dialog(tab="Advanced", group="Conservation of mass, energy, momentum"));
  parameter Boolean kineticTerm 
      " = true, include kinetic term in momentum balance" 
                                annotation(Evaluate=true, Dialog(tab="Advanced", group="Conservation of mass, energy, momentum"));
  parameter Real K1=1 
      "Enhancement factor for heat transfer area pipe 1(=>parallel tubes)" 
                                                                          annotation(Dialog(tab="General", group="Heat transfer"));
  parameter Real K2=1 
      "Enhancement factor for heat transfer area pipe 2(=>parallel tubes)" 
                                                                          annotation(Dialog(tab="General", group="Heat transfer"));
  parameter Boolean dynamicTerm=false 
      " = true, include dynamic term in momentum balance, only if not lumped_dp and not static"
                                                                                              annotation(Evaluate=true, Dialog(tab="Advanced", group="Conservation of mass, energy, momentum"));
    
  //Pressure drop and heat transfer    
  replaceable package WallFriction = 
      PressureLosses.Utilities.WallFriction.QuadraticTurbulent extends 
      Modelica_Fluid.PressureLosses.Utilities.WallFriction.PartialWallFriction 
      "Characteristic of wall friction"                                                          annotation(choicesAllMatching, Dialog(tab="General", group="Pressure drop"));
  parameter SI.Length roughness_1=2.5e-5 
      "Absolute roughness of pipe (default = smooth steel pipe)" 
                                                               annotation(Dialog(tab="General", group="Pressure drop"));
  parameter SI.Length roughness_2=2.5e-5 
      "Absolute roughness of pipe (default = smooth steel pipe)" 
                                                               annotation(Dialog(tab="General", group="Pressure drop"));
  parameter SI.DynamicViscosity eta_nominal_M1=0.01 
      "Nominal dynamic viscosity of medium 1(e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
                                                                                           annotation(Dialog(tab="General", group="Pressure drop"));
  parameter SI.DynamicViscosity eta_nominal_M2=0.01 
      "Nominal dynamic viscosity of medium 1(e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
                                                                                       annotation(Dialog(tab="General", group="Pressure drop"));
  parameter Boolean use_eta_nominal=false 
      "= true, if eta_nominal is used, otherwise computed from medium" 
                                                                     annotation(Evaluate=true, Dialog(tab="General", group="Pressure drop"));
  replaceable model HeatTransfer_1 = 
      Modelica_Fluid.HeatTransfer.PipeHT_constAlpha 
      extends Modelica_Fluid.HeatTransfer.PartialPipeHeatTransfer 
      "Heat transfer model"                                                                       annotation(choicesAllMatching, Dialog(tab="General", group="Heat transfer"));
  replaceable model HeatTransfer_2 = 
      Modelica_Fluid.HeatTransfer.PipeHT_constAlpha 
      extends Modelica_Fluid.HeatTransfer.PartialPipeHeatTransfer 
      "Heat transfer model"                                                                       annotation(choicesAllMatching, Dialog(tab="General", group="Heat transfer"));
  //Display variables
  SI.HeatFlowRate Q_flow_1 "Total heat flow rate of inner pipe";
  SI.HeatFlowRate Q_flow_2 "Total heat flow rate of outer pipe";
    
  Modelica_Fluid.Components.LongPipe pipe_1(
    redeclare package Medium = Medium_1,
    n=n,
    static=static,
    lumped_dp=lumped_dp,
    kineticTerm=kineticTerm,
    crossSectionType=Modelica_Fluid.Types.CrossSectionTypes.Circular,
    length=length,
    use_wall=true,
    area_h=Modelica.Constants.pi*di_1*length*K1,
    redeclare HeatTransfer_1 heat,
    initOption=initOption_1,
    use_T_start=use_T_start_1,
    T_start=T_start_1,
    h_start=h_start_1,
    X_start=X_start_1,
    mflow_start=mflow_start_1,
    d_inner=di_1,
    d_outer=da_1,
    redeclare package WallFriction = WallFriction,
    roughness=roughness_1,
    use_eta_nominal=use_eta_nominal,
    eta_nominal=eta_nominal_M1,
    redeclare final model Wall = 
        Modelica_Fluid.Components.Wall_constProps (
        d_wall=d_wall,
        c_wall=c_wall,
        initOption=if initWall_steadyState then 3 else 2,
        T_start=T_start_wall),
      dynamicTerm=dynamicTerm) 
                             annotation (extent=[-40,-60; 20,0]);
    
  Modelica_Fluid.Components.LongPipe pipe_2(
    redeclare package Medium = Medium_2,
    n=n,
    static=static,
    lumped_dp=lumped_dp,
    kineticTerm=kineticTerm,
    crossSectionType=Modelica_Fluid.Types.CrossSectionTypes.General,
    length=length,
    redeclare HeatTransfer_2 heat,
    use_T_start=use_T_start_2,
    T_start=T_start_2,
    h_start=h_start_2,
    X_start=X_start_2,
    initOption=initOption_2,
    mflow_start=mflow_start_2,
    P_inner=Modelica.Constants.pi*(da_1 + di_2),
    A_inner=Modelica.Constants.pi/4*(di_2*di_2 - da_1*da_1),
    area_h=Modelica.Constants.pi*da_1*length*K2,
    p_a_start=p_a_start1,
    p_b_start=p_b_start2,
    redeclare package WallFriction = WallFriction,
    roughness=roughness_2,
    use_eta_nominal=use_eta_nominal,
    eta_nominal=eta_nominal_M2,
    dynamicTerm=dynamicTerm) 
              annotation (extent=[-40,88; 20,28]);
  annotation (Diagram(Line(points=[-10,36; -10,-8], style(
          color=1,
          rgbcolor={255,0,0},
          thickness=2))), Icon(
      Rectangle(extent=[-100,-26; 100,-30], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95},
          fillPattern=7)),
      Rectangle(extent=[-100,30; 100,26], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95},
          fillPattern=7)),
      Rectangle(extent=[-100,60; 100,30], style(
          color=0,
          rgbcolor={0,0,0},
          gradient=2,
          fillColor=70,
          rgbfillColor={0,63,125})),
      Rectangle(extent=[-100,-30; 100,-60], style(
          color=0,
          rgbcolor={0,0,0},
          gradient=2,
          fillColor=70,
          rgbfillColor={0,63,125})),
      Rectangle(extent=[-100,26; 100,-26], style(
          color=69,
          rgbcolor={0,128,255},
          gradient=2,
          fillColor=69,
          rgbfillColor={0,128,255})),
      Text(
        extent=[-100,-60; 100,-100],
        string="%name",
        style(color=3, rgbcolor={0,0,255}))));
  Modelica_Fluid.Interfaces.FluidPort_b port_b1(redeclare package Medium = 
        Medium_1) annotation (extent=[100,-12; 120,8]);
  Modelica_Fluid.Interfaces.FluidPort_a port_a1(redeclare package Medium = 
        Medium_1) annotation (extent=[-120,-12; -100,8]);
  Modelica_Fluid.Interfaces.FluidPort_a port_a2(redeclare package Medium = 
        Medium_2) annotation (extent=[-120,36; -100,56]);
  Modelica_Fluid.Interfaces.FluidPort_b port_b2(redeclare package Medium = 
        Medium_2) annotation (extent=[100,-56; 120,-36]);
    
equation 
  Q_flow_1=sum(pipe_1.Qs_flow);
  Q_flow_2=sum(pipe_2.Qs_flow);
  connect(pipe_2.thermalPort, pipe_1.thermalPort);
  connect(pipe_2.port_b, port_b2) annotation (points=[20,58; 60,58; 60,-46; 110,
        -46], style(
      color=69,
      rgbcolor={0,127,255},
      thickness=2,
      gradient=2,
      fillColor=42,
      rgbfillColor={213,0,0}));
  connect(pipe_1.port_b, port_b1) annotation (points=[20,-30; 42,-30; 42,-2;
        110,-2], style(
      color=69,
      rgbcolor={0,127,255},
      thickness=2,
      gradient=2,
      fillColor=42,
      rgbfillColor={213,0,0}));
  connect(pipe_1.port_a, port_a1) annotation (points=[-40.6,-30; -75.3,-30; 
          -75.3,-2; -110,-2],
                            style(
      color=69,
      rgbcolor={0,127,255},
      thickness=2,
      gradient=2,
      fillColor=42,
      rgbfillColor={213,0,0}));
  connect(pipe_2.port_a, port_a2) annotation (points=[-40.6,58; -76,58; -76,46; 
          -110,46],
                  style(
      color=69,
      rgbcolor={0,127,255},
      thickness=2,
      gradient=2,
      fillColor=42,
      rgbfillColor={213,0,0}));
end HeatExchanger;
  
end Components;
