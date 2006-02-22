package FluidStorage 
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
        package Medium=Medium, m_flow(each start=0), mXi_flow(each start=0)) 
    annotation (extent=[-30,100; 30,108]);
  Modelica_Fluid.Interfaces.FluidPort_ArrayIcon bottomPorts[n_bottomPorts](redeclare 
        package Medium=Medium, m_flow(each start=0), mXi_flow(each start=0)) 
    annotation (extent=[-30,-108; 30,-100],   rotation=90);
  Modelica_Fluid.Interfaces.FluidPort_ArrayIcon sidePorts[n_sidePorts](redeclare 
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
</HTML>", revisions="<html>
<ul>
<li><i>1 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from a previous version of Modelica_Fluid</li>
</ul>
</html>"),
    Diagram);
end Tank;
end FluidStorage;
