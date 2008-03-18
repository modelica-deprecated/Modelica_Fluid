within Modelica_Fluid.WorkInProgress;
package Components 
  
model IsolatedPipe 
    "Model of an isolated pipe consisting of n pipe segments/FiniteVolumes" 
    
  replaceable package Medium = 
    Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
      annotation (choicesAllMatching = true);
    
  extends Modelica_Fluid.Interfaces.PartialInitializationParameters;
    
  parameter Integer nVolumes(min=1)=1 "Number of pipe segments/finite volumes";
    
  parameter SI.Length L "Length of pipe";
  parameter SI.AbsolutePressure dp_nominal(min=1.e-10) = 1 
      "|frictionType = ConstantLaminar or ConstantTurbulent| Nominal pressure drop";
    
  parameter SI.MassFlowRate m_flow_nominal = 1E-3 
      "Nominal mass flow rate at nominal pressure drop";
    
  parameter SI.Area A_a;
  parameter SI.Area A_b=A_a;
    
  parameter SI.Length Z_a=0;
  parameter SI.Length Z_b=Z_a;
    
  parameter Boolean dynamicMomentumBalance=false 
      "If false, der(m_flow) is neglected in momentum balance" 
                                                 annotation(Evaluate=true,
      Dialog(tab="Level of Detail"));
  parameter Boolean includeKineticTerm=false 
      "If false, d*v^2 is neglected in momentum balance" 
                                             annotation(Evaluate=true,
      Dialog(tab="Level of Detail"));
  parameter Boolean includeViscosity=false 
      "If false, artifical viscosity is neglected" 
                                          annotation(Evaluate=true, Dialog(tab=
          "Level of Detail"));
  parameter Real viscosityFactor1=0 annotation(Dialog(enable=includeViscosity,tab="Level of Detail"));
  parameter Real viscosityFactor2=1 annotation(Dialog(enable=includeViscosity,tab="Level of Detail"));
    
  Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare model Medium = Medium) 
              annotation (extent=[-120, -10; -100, 10]);
  Modelica_Fluid.Interfaces.FluidPort_b port_b(redeclare model Medium = Medium) 
              annotation (extent=[120, -10; 100, 10]);
  Modelica_Fluid.WorkInProgress.Utilities.PipeSegment pipeSegment[nVolumes](
      redeclare package Medium = Medium,
      each initType = initType,
      each p_start = p_start,
      each use_T_start = use_T_start,
      each T_start = T_start,
      each h_start = h_start,
      each X_start = X_start,
      each L=L/nVolumes,
      each dp_nominal=dp_nominal/nVolumes,
      each A_a=A_a "has to be corrected: linear distribution of A",
      each Z_a=Z_a "has to be corrected: linear distribution of Z",
      each m_flow_nominal=m_flow_nominal,
      each dynamicMomentumBalance=dynamicMomentumBalance,
      each includeKineticTerm=includeKineticTerm,
      each includeViscosity=includeViscosity,
      each viscosityFactor1=viscosityFactor1,
      each viscosityFactor2=viscosityFactor2);
    
annotation (Icon(
    Rectangle(extent=[-100, 60; 100, -60], style(color=0, fillColor=8)),
    Rectangle(extent=[-100, 34; 100, -36], style(
        color=69,
        gradient=2,
        fillColor=69)),
    Text(
      extent=[-150, 125; 150, 65],
      string="%name",
      style(gradient=2, fillColor=69)),
      Ellipse(extent=[-58,14; -28,-14], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0})),
      Ellipse(extent=[22,14; 52,-14], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0}))));
equation 
  connect(port_a, pipeSegment[1].port_a);
  connect(port_b, pipeSegment[nVolumes].port_b);
  for i in 1:nVolumes - 1 loop
    connect(pipeSegment[i].port_b, pipeSegment[i + 1].port_a);
  end for;
end IsolatedPipe;
  
  model ShortPipe2 
    "Short pipe with two volumes, wall friction and gravity effect" 
    
    extends Modelica_Fluid.Interfaces.PartialTwoPortTransport;
    replaceable package WallFriction = 
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent
      extends 
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.PartialWallFriction
      "Characteristic of wall friction"  annotation(choicesAllMatching=true);
    
    parameter SI.Length length "Length of pipe";
    parameter SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
    parameter SI.Length height_ab = 0.0 "Height of port_b over port_a" annotation(Evaluate=true);
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
  /*
  parameter Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.Temp 
    flowDirection=
            Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.UseGlobalFluidOption 
    "Unidirectional (port_a -> port_b) or bidirectional flow component" 
     annotation(Dialog(tab="Advanced"));
*/
    parameter SI.AbsolutePressure dp_small = 1 
      "Turbulent flow for wall friction if |dp| >= dp_small" 
      annotation(Dialog(tab="Advanced", enable=WallFriction.use_dp_small));
    final parameter SI.Volume V = Modelica.Constants.pi*(diameter/2)^2*length;
    
    parameter Modelica_Fluid.Types.Init.Temp initVolume1=
               Modelica_Fluid.Types.Init.NoInit 
      "Initialization option for volume 1" 
      annotation(Dialog(tab = "Initialization"));
    
    parameter Modelica_Fluid.Types.Init.Temp initVolume2=
              Modelica_Fluid.Types.Init.NoInit 
      "Initialization option for volume 2" 
      annotation(Dialog(tab = "Initialization"));
    
  /*
  parameter Medium.AbsolutePressure p_start = Medium.p_default 
    "Start value of pressure" 
    annotation(Dialog(tab = "Initialization"));
  parameter Boolean use_T_start = true "Use T_start if true, otherwise h_start"
    annotation(Dialog(tab = "Initialization"), Evaluate=true);
  parameter Medium.Temperature T_start=
    if use_T_start then Medium.T_default else Medium.temperature_phX(p_start,h_start,X_start) 
    "Start value of temperature" 
    annotation(Dialog(tab = "Initialization", enable = use_T_start));
  parameter Medium.SpecificEnthalpy h_start=
    if use_T_start then Medium.specificEnthalpy_pTX(p_start, T_start, X_start) else Medium.h_default 
    "Start value of specific enthalpy" 
    annotation(Dialog(tab = "Initialization", enable = not use_T_start));
  parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
    "Start value of mass fractions m_i/m" 
    annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
*/
    
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
        Ellipse(extent=[-90,15; -60,-15], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=0,
              rgbfillColor={0,0,0})),
        Ellipse(extent=[60,15; 90,-15],   style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=0,
              rgbfillColor={0,0,0}))),       Documentation(info="<html>
<p>
Simple pipe model consisting of two volumes, 
wall friction (with different friction correlations)
and gravity effect. This model is mostly used to demonstrate how
to build up more detailed models from the basic components.
</p>
</html>"),
      Diagram,
      Coordsys(grid=[1,1], scale=0));
    PressureLosses.WallFrictionAndGravity frictionAndGravity(
      redeclare package Medium = Medium,
      flowDirection=flowDirection,
      redeclare package WallFriction = WallFriction,
      diameter=diameter,
      roughness=roughness,
      use_nominal=use_nominal,
      eta_nominal=eta_nominal,
      d_nominal=d_nominal,
      from_dp=true,
      dp_small=dp_small,
      show_Re=false,
      length=length,
      height_ab=height_ab) 
                         annotation (extent=[-10,-10; 10,10]);
    Modelica_Fluid.Pipes.BaseClasses.PortVolume volume1(
      redeclare package Medium = Medium,
      p_start=p_a_start,
      use_T_start=use_T_start,
      T_start=T_start,
      h_start=h_start,
      X_start=X_start,
      V=V/2,
      initType=initVolume1) 
      annotation (extent=[-70,-10; -50,10]);
    Modelica_Fluid.Pipes.BaseClasses.PortVolume volume2(
      redeclare package Medium = Medium,
      p_start=p_b_start,
      use_T_start=use_T_start,
      T_start=T_start,
      h_start=h_start,
      X_start=X_start,
      V=V/2,
      initType=initVolume2) 
      annotation (extent=[50,-10; 70,10]);
  equation 
    connect(volume1.port, port_a) 
      annotation (points=[-60,0; -100,0], style(color=69, rgbcolor={0,127,255}));
    connect(volume1.port, frictionAndGravity.port_a) 
      annotation (points=[-60,0; -10,0], style(color=69, rgbcolor={0,127,255}));
    connect(frictionAndGravity.port_b, volume2.port) 
      annotation (points=[10,0; 60,0], style(color=69, rgbcolor={0,127,255}));
    connect(volume2.port, port_b) 
      annotation (points=[60,0; 100,0], style(color=69, rgbcolor={0,127,255}));
  end ShortPipe2;
  
model OpenTank "Tank with three inlet/outlet-arrays at variable heights" 
      import Modelica_Fluid.Types;
    import Modelica_Fluid;
  replaceable package Medium = 
    Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
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
    
  parameter Medium.AbsolutePressure p_ambient=ambient.default_p_ambient 
      "Tank surface pressure" 
    annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));
  parameter Medium.Temperature T_ambient = ambient.default_T_ambient 
      "Tank surface Temperature" 
    annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));
    
  parameter Types.Init.Temp initType=
            Types.Init.InitialValues "Initialization option" 
    annotation(Evaluate=true,Dialog(tab = "Ambient and Initialization", group = "Initialization"));
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
    
  Modelica_Fluid.WorkInProgress.FluidStorage.FluidPort_ArrayIcon topPorts[
                                                                       n_topPorts](redeclare 
        package Medium=Medium, m_flow(each start=0), mXi_flow(each start=0)) 
    annotation (extent=[-30,100; 30,108]);
  Modelica_Fluid.WorkInProgress.FluidStorage.FluidPort_ArrayIcon bottomPorts[
                                                                          n_bottomPorts](redeclare 
        package Medium=Medium, m_flow(each start=0), mXi_flow(each start=0)) 
    annotation (extent=[-30,-108; 30,-100],   rotation=90);
  Modelica_Fluid.WorkInProgress.FluidStorage.FluidPort_ArrayIcon sidePorts[
                                                                        n_sidePorts](redeclare 
        package Medium=Medium, m_flow(each start=0), mXi_flow(each start=0)) 
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
    
Modelica_Fluid.Ambient ambient;
  protected 
  Real Q_lost = 0 "Wärmeverlust (currently dummy)";
  parameter Medium.SpecificEnthalpy h_start = Medium.specificEnthalpy_pTX(p_ambient, T_start, X_start);
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
  Modelica_Fluid.WorkInProgress.FluidStorage.TankAttachment tankAttachmentTop[
                                                                           n_topPorts](
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
  Modelica_Fluid.WorkInProgress.FluidStorage.TankAttachment 
      tankAttachmentBottom[                                                   n_bottomPorts](
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
  Modelica_Fluid.WorkInProgress.FluidStorage.TankAttachment tankAttachmentSide[
                                                                            n_sidePorts](
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
  if initType == Types.Init.NoInit then
    // no initial equations
  elseif initType == Types.Init.InitialValues then
    level = level_start;
    medium.T = T_start;
    medium.Xi = X_start[1:Medium.nXi];
  elseif initType == Types.Init.SteadyState then
    der(level) = 0;
    der(medium.h) = 0;
    der(medium.Xi) = zeros(Medium.nXi);
  elseif initType == Types.Init.SteadyStateHydraulic then
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
 
</HTML>",
        revisions="<html>
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
  
model Tank "Obsolet component (use instead Components.OpenTank)" 
    import Modelica_Fluid.Types;
    import Modelica_Fluid;
  replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
    annotation (choicesAllMatching=true);
  parameter SI.Area area "Tank area";
  parameter SI.Area pipeArea "Area of outlet pipe";
  parameter SI.Volume V0 = 0 "Volume of the liquid when the level is zero";
  parameter SI.Height H0 = 0 
      "Height of zero level reference over the bottom port";
  parameter Medium.AbsolutePressure p_ambient=ambient.default_p_ambient 
      "Tank surface pressure";
  parameter Types.Init.Temp initType=
            Types.Init.NoInit "Initialization option" 
    annotation(Dialog(tab = "Initialization"));
  parameter Boolean use_T_start = true "Use T_start if true, otherwise h_start"
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
    
  Modelica_Fluid.Interfaces.FluidPort_b port(
      redeclare package Medium = Medium,
      m_flow(start=0),
      mXi_flow(each start=0)) 
    annotation (extent=[-10, -120; 10, -100], rotation=90);
  Medium.BaseProperties medium(
    preferredMediumStates=true,
    p(start=p_ambient),
    T(start=T_start),
    Xi(start=X_start[1:Medium.nXi]));
    outer Modelica_Fluid.Ambient ambient "Ambient conditions";
    
  SI.Height level(start=level_start,stateSelect=StateSelect.prefer) 
      "Height of tank level over the zero reference";
  Medium.AbsolutePressure p_bottom "Pressure at bottom of tank";
  SI.Energy U "Internal energy of tank volume";
  SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";
  Real m(quantity=Medium.mediumName, unit="kg", stateSelect=StateSelect.never) 
      "Mass of tank volume";
  Real mX[Medium.nX](quantity=Medium.substanceNames, each unit="kg") 
      "Component masses of the independent substances";
    
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
  p_bottom = (medium.d*ambient.g*(level+H0)) + p_ambient;
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
  if initType == Types.Init.NoInit then
    // no initial equations
  elseif initType == Types.Init.InitialValues then
    level = level_start;
    if use_T_start then
      medium.T = T_start;
    else
      medium.h = h_start;
    end if;
    medium.Xi = X_start[1:Medium.nXi];
  elseif initType == Types.Init.SteadyState then
    der(level) = 0;
    der(medium.h) = 0;
    der(medium.Xi) = zeros(Medium.nXi);
  elseif initType == Types.Init.SteadyStateHydraulic then
    der(level) = 0;
    if use_T_start then
      medium.T = T_start;
    else
      medium.h = h_start;
    end if;
    medium.Xi = X_start[1:Medium.nXi];
  else
    assert(false,"Unsupported initialization option initType = " + String(initType)
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
</HTML>", revisions="<html>
<ul>
<li><i>1 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from a previous version of Modelica_Fluid</li>
</ul>
</html>"),
    Diagram);
end Tank;
  
    model PortVolume 
    "Fixed volume associated with a port by the finite volume method (used to build up physical components; fulfills mass and energy balance)" 
    import SI = Modelica.SIunits;
    import Modelica_Fluid.Types;
      extends 
      Modelica_Fluid.WorkInProgress.Interfaces.PartialInitializationParameters;
    
      replaceable package Medium = 
          Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
          annotation (choicesAllMatching = true);
    
      parameter SI.Volume V "Volume";
    
      Modelica_Fluid.Interfaces.FluidPort_a port(
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
         zeros(Medium.nC) = port.mC_flow "Trace substances";
    
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

  model LumpedPipe 
    "Short pipe with one volume, wall friction and gravity effect" 
    replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
      annotation (choicesAllMatching = true);
    parameter Modelica_Fluid.Types.Init.Temp initType=Modelica_Fluid.Types.Init.NoInit 
      "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_a_start 
      "Start value of pressure at port_a" 
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_b_start = p_a_start 
      "Start value of pressure at port_b" 
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
      extends 
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
    parameter Modelica_Fluid.Types.FlowDirection.Temp flowDirection=
        Modelica_Fluid.Types.FlowDirection.Bidirectional 
      "Unidirectional (port_a -> port_b) or bidirectional flow component" 
       annotation(Dialog(tab="Advanced"));
    parameter Modelica.SIunits.AbsolutePressure dp_small=1 
      "Turbulent flow if |dp| >= dp_small (only used if WallFriction=QuadraticTurbulent)"
      annotation(Dialog(tab="Advanced", enable=WallFriction.use_dp_small));
    
      Modelica_Fluid.Interfaces.FluidPort_a port_a(
                                    redeclare package Medium = Medium) 
      "Fluid connector a (positive design flow direction is from port_a to port_b)"
        annotation (extent=[-110,-10; -90,10]);
      Modelica_Fluid.Interfaces.FluidPort_b port_b(
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
      Coordsys(grid=[1,1], scale=0),
      uses(Modelica_Fluid(version="1.0 Beta 2"), Modelica(version="2.2.2")));
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
        show_Re=false,
        p_a_start=p_a_start,
        p_b_start=(p_a_start+p_b_start)/2,
        T_start=T_start,
        h_start=h_start,
        X_start=X_start) annotation (extent=[-60,-10; -40,10]);
    Modelica_Fluid.WorkInProgress.Components.PortVolume volume(
      redeclare package Medium = Medium,
      V=Modelica.Constants.pi*(diameter/2)^2*length,
      initType=initType,
      p_start=(p_a_start+p_b_start)/2,
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
        show_Re=false,
        p_a_start=(p_a_start+p_b_start)/2,
        p_b_start=p_b_start,
        T_start=T_start,
        h_start=h_start,
        X_start=X_start,
        use_T_start=true) 
                         annotation (extent=[40,-10; 60,10]);
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
end Components;
