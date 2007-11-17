package Volumes "Generic volume, tank and other volume type components" 
    model PortVolume "Volume in a connection point" 
    
    annotation (
      Icon(
        Ellipse(extent=[-100,100; 100,-100], style(
            color=0,
            rgbcolor={0,0,0},
            gradient=3,
            fillColor=68,
            rgbfillColor={170,213,255})),
        Text(extent=[-144,178; 146,116], string="%name"),
        Text(
          extent=[-130,-108; 144,-150],
          style(color=0),
          string="V=%V")),
      Documentation(info="<html>
Ideally mixed volume of constant size with two fluid ports and one medium model. The flow properties are computed from the upstream quantities, pressures are equal in both nodes and the medium model. Heat transfer through a thermal port is possible, it equals zero if the port remains unconnected. The thermal port temperature is equal to the medium temperature.
</html>"),
      Diagram);
    
      replaceable package Medium = 
        Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
          annotation (choicesAllMatching = true);
      parameter SI.Volume V "Volume";
      parameter Types.Init.Temp initType=
                Types.Init.NoInit "Initialization option" 
        annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p_start = Medium.p_default 
      "Start value of pressure" 
        annotation(Dialog(tab = "Initialization"));
      parameter Boolean use_T_start = true 
      "= true, use T_start, otherwise h_start" 
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
    
      Interfaces.FluidStatePort_a port(redeclare package Medium = Medium) 
      "Port inside the volume" 
                         annotation (extent=[-10,-10; 10,10]);
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a thermalPort 
      "Thermal port" 
      annotation (extent=[-20,88; 20,108]);
    
      Medium.BaseProperties medium(
        preferredMediumStates=true,
        p(start=p_start),
        h(start=h_start),
        T(start=T_start),
        Xi(start=X_start[1:Medium.nXi]));
      SI.Energy U "Internal energy of fluid";
      SI.Mass m "Mass of fluid";
      SI.Mass mXi[Medium.nXi] "Masses of independent components in the fluid";
    
    equation 
     thermalPort.T = medium.T;
    
    // boundary conditions
      port.p  = medium.p;
      port.h  = medium.h;
      port.Xi = medium.Xi;
    
    // Total quantities
      m   = V*medium.d;
      mXi = m*medium.Xi;
      U   = m*medium.u;
    
    // Mass and energy balance
      der(m) = port.m_flow;
      der(mXi) = port.mXi_flow;
      der(U) = port.H_flow + thermalPort.Q_flow;
      zeros(Medium.nC) = port.mC_flow;
    
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
   extends Modelica_Fluid.Icons.VariantLibrary;
  
    model MixingVolume 
    "Mixing volume with inlet and outlet ports (flow reversal is allowed)" 
    extends Modelica_Fluid.Volumes.BaseClasses.PartialLumpedVolume(
                                                          V_lumped=V, Ws_flow=0);
    parameter SI.Volume V "Volume";
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a thermalPort 
      "Thermal port" 
      annotation (extent=[-20,88; 20,108]);
    annotation (
      Icon(
        Ellipse(extent=[-100,100; 100,-100], style(
            color=0,
            rgbcolor={0,0,0},
            gradient=3,
            fillColor=68,
            rgbfillColor={170,213,255})),
        Text(extent=[-144,178; 146,116], string="%name"),
        Text(
          extent=[-130,-108; 144,-150],
          style(color=0),
          string="V=%V")),
      Documentation(info="<html>
Ideally mixed volume of constant size with two fluid ports and one medium model. The flow properties are computed from the upstream quantities, pressures are equal in both nodes and the medium model. Heat transfer through a thermal port is possible, it equals zero if the port remains unconnected. The thermal port temperature is equal to the medium temperature.
</html>"),
      Diagram);
    equation 
    thermalPort.T = medium.T;
    Qs_flow = thermalPort.Q_flow;
    end MixingVolume;
  
  model SweptVolume 
    "varying cylindric volume depending on the postition of the piston" 
    extends BaseClasses.PartialLumpedVolume;
    parameter SI.Area pistonCrossArea "cross sectional area of pistion";
    parameter SI.Volume clearance "remaining volume at zero piston stroke";
    annotation (Diagram, Icon(
        Rectangle(extent=[-44,62; 44,-30], style(
            pattern=0,
            thickness=4,
            fillColor=31,
            rgbfillColor={170,213,255},
            fillPattern=1)),
        Polygon(points=[-40,60; -40,44; -44,44; -44,64; 44,64; 44,44; 40,44; 40,
              60; -40,60], style(
            color=10,
            rgbcolor={95,95,95},
            smooth=0,
            fillColor=9,
            rgbfillColor={135,135,135},
            fillPattern=8)),
        Polygon(points=[-44,34; -40,34; -40,-60; -44,-60; -44,34], style(
            color=10,
            rgbcolor={95,95,95},
            smooth=0,
            fillColor=9,
            rgbfillColor={135,135,135},
            fillPattern=8)),
        Polygon(points=[40,34; 44,34; 44,-60; 40,-60; 40,34], style(
            color=10,
            rgbcolor={95,95,95},
            smooth=0,
            fillColor=9,
            rgbfillColor={135,135,135},
            fillPattern=8)),
        Rectangle(extent=[-40,-30; 40,-40], style(
            color=10,
            rgbcolor={95,95,95},
            fillColor=9,
            rgbfillColor={135,135,135},
            fillPattern=7)),
        Rectangle(extent=[-6,-40; 6,-92], style(
            color=10,
            rgbcolor={95,95,95},
            fillColor=9,
            rgbfillColor={135,135,135},
            fillPattern=7)),
        Line(points=[-102,0; -70,0; -70,40; -44,40], style(
            color=31,
            rgbcolor={170,213,255},
            smooth=0)),
        Line(points=[44,40; 70,40; 70,0; 100,0], style(
            color=31,
            rgbcolor={170,213,255},
            smooth=0))),
      Documentation(info="<html>
<p> Mixing volume with varying size. The size of the volume is given by:</p>
<ul>
  <li>cross sectional piston area</li>
  <li>piston stroke given by the flange position s</li>
  <li>clearance (volume at flang position = 0)</li>
</ul> 
 
<p> The flange position has to be equal or greater than zero. Otherwise the simulation stops. The force of the flange results from the pressure difference between medium and ambient pressure and the cross sectional piston area. For using the component, a top level instance of the ambient model with the inner attribute is needed.</p>
<p> The pressure at both fluid ports equals the medium pressure in the volume. No suction nor discharge valve is included in the model.</p>
<p>The thermal port is directly connected to the medium. The temperature of the thermal port equals the medium temperature. The heat capacity of the cylinder and the piston are not includes in the model.</p>
</html>",
        revisions="<html>
<ul>
<li><i>29 Oct 2007</i>
    by Carsten Heinrich:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
    Modelica.Mechanics.Translational.Interfaces.Flange_b flange 
      "translation flange for piston" annotation (extent=[-10,-110; 10,-90]);
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a thermalPort 
      "Thermal port" 
    annotation (extent=[-20,88; 20,108]);
    outer Modelica_Fluid.Ambient ambient "Ambient conditions";
    
  equation 
    assert(flange.s >= 0, "Piston stroke (given by flange.s) must not be smaller than zero!");
    
    // volume size
    V_lumped = clearance + flange.s * pistonCrossArea;
    
    flange.f = (medium.p - ambient.default_p_ambient) * pistonCrossArea;
    thermalPort.T = medium.T;
    
    // energy balances
    Ws_flow = medium.p * pistonCrossArea * (-der(flange.s));
    Qs_flow = thermalPort.Q_flow;
  end SweptVolume;
  
model OpenTank "Open tank with inlet/outlet ports at the bottom" 
    replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
    annotation (choicesAllMatching=true);
    
//Transformation of kinetic energy
    parameter Boolean p_static_at_port=false 
      "=true, kinetic energy and dissipation is accounted for in port pressure"
                                                                                                        annotation(Evaluate=true, Dialog(tab="Advanced"));
    parameter Real[n_ports] zeta_in=fill(0, n_ports) 
      "Hydraulic resistance into tank, 1 for total dissipation of kinetic energy and uniform flow distribution in pipe"
                                                                                                        annotation(Dialog(tab="Advanced",enable=p_static_at_pot==false));
    parameter Real[n_ports] zeta_out=fill(1, n_ports) 
      "Hydraulic resistance out of tank, 0 for ideal smooth outlet" 
                                                                  annotation(Dialog(tab="Advanced",enable=p_static_at_pot==false));
    
//Tank geometry  
    parameter SI.Height height "Height of tank";
    parameter SI.Area area "Area of tank";
    parameter SI.Volume V0=0 "Volume of the liquid when the level is zero";
    
//Port definitions 
    parameter Integer n_ports(min=1) = 1 "Number of bottom ports (min=1)" 
     annotation(Dialog(group="bottomPorts (= pipes at bottom of tank; in and out flow of tank)"));
    parameter SI.Diameter pipe_diameters[n_ports] 
      "Inner (hydraulic) diameters of bottom ports (array)" 
     annotation(Dialog(group="bottomPorts (= pipes at bottom of tank; in and out flow of tank)", enable=n_bottomPorts > 0));
    Interfaces.FluidPorts_b ports[n_ports](
    redeclare package Medium = Medium,
    m_flow(each start=0),
    mXi_flow(each start=0)) 
    annotation (extent=[-5,-125; 5,-85], rotation=90);
    
//Ambient  
    parameter Medium.AbsolutePressure p_ambient=ambient.default_p_ambient 
      "Tank surface pressure" 
    annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));
    parameter Medium.Temperature T_ambient=ambient.default_T_ambient 
      "Tank surface Temperature" 
    annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));
    outer Modelica_Fluid.Ambient ambient;
    
//Initialization
    parameter Types.Init.Temp initType=Types.Init.InitialValues 
      "Initialization option" 
    annotation(Evaluate=true,Dialog(tab = "Ambient and Initialization", group = "Initialization"));
    parameter SI.Height level_start "Start value of tank level" 
    annotation(Dialog(tab="Ambient and Initialization", group = "Initialization"));
    parameter Medium.Temperature T_start=T_ambient "Start value of temperature"
    annotation(Dialog(tab = "Ambient and Initialization", group = "Initialization"));
    parameter Medium.MassFraction X_start[Medium.nX]=Medium.X_default 
      "Start value of mass fractions m_i/m" 
    annotation (Dialog(tab="Ambient and Initialization", group = "Initialization", enable=Medium.nXi > 0));
    
//Tank properties  
    Medium.BaseProperties medium(
      preferredMediumStates=true,
      p(start=p_ambient),
      T(start=T_start),
      h(start=h_start),
      Xi(start=X_start[1:Medium.nXi]));
    SI.Height level(stateSelect=StateSelect.default, start=level_start) 
      "Level height of tank";
    SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";
    SI.Energy U "Internal energy of tank volume";
    SI.Mass m "Mass of fluid in tank";
    SI.Mass mXi[Medium.nXi] "Masses of independent components in the fluid";
    SI.Pressure p_static "bottom tank pressure";
    
  protected 
    parameter Medium.SpecificEnthalpy h_start=Medium.specificEnthalpy_pTX(
        p_ambient,
        T_start,
        X_start);
    parameter SI.Area[n_ports] pipeArea=Modelica.Constants.pi/4*{pipe_diameters[
        i]^2 for i in 1:n_ports};
    
equation 
  //Total quantities
    V = area*level + V0 "Volume of fluid";
    m = V*medium.d "Mass of fluid";
    mXi = m*medium.Xi "Mass of fluid components";
    U = m*medium.u "Internal energy of fluid";
    medium.p = p_ambient;
    
  // Mass balances
    der(m) = sum(ports.m_flow);
    for i in 1:Medium.nXi loop
      der(mXi[i]) = sum(ports[:].mXi_flow[i]);
    end for;
    
  // Energy balance
    if Medium.singleState then
      der(U) = sum(ports.H_flow);
                               //Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)
    else
      der(U) = sum(ports.H_flow) - p_ambient*der(V);
    end if;
    assert(level <= height, "Tank is full (level = height = " + String(level) + ")");
    assert(level > 0, "Tank is empty (level = 0), tank model is not designed to allow air flow through ports");
    
//Determine port properties  
    p_static = level*ambient.g*medium.d + p_ambient;
    for i in 1:n_ports loop
      ports[i].H_flow = semiLinear(
        ports[i].m_flow,
        ports[i].h,
        medium.h);
      ports[i].mXi_flow = semiLinear(
        ports[i].m_flow,
        ports[i].Xi,
        medium.Xi);
      if p_static_at_port then
        ports[i].p = p_static;
      else
       ports[i].p = p_static - smooth(2, noEvent(ports[i].m_flow^2/(2*medium.d*
          pipeArea[i]^2)*(if ports[i].m_flow < 0 then (1 + zeta_out[i]) else (1
           - zeta_in[i]))));
       end if;
    end for;
    
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
    annotation (
      Icon(
        Rectangle(extent=[-100,100; 100,0], style(color=7, fillColor=7)),
        Rectangle(extent=DynamicSelect([-100,-100; 100,10], [-100,-100; 100,(-100
               + 200*level/height)]), style(
            color=69,
            rgbcolor={0,127,255},
            fillColor=71,
            rgbfillColor={85,170,255},
            fillPattern=1)),
        Line(points=[-100,100; -100,-100; 100,-100; 100,100], style(
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
          extent=[-129,53; 130,39],
          style(color=0),
          Line(points=[-100,100; 100,100], style(
              color=0,
              rgbcolor={0,0,0},
              pattern=3)),
        string="start = %level_start m"),
        Text(
          extent=[-95,30; 95,5],
          style(color=0),
          string=DynamicSelect(" ", realString(
              level,
              1,
              integer(precision)))),
        Line(points=[-100,100; 100,100], style(
            color=0,
            rgbcolor={0,0,0},
            pattern=3))),
      Documentation(info="<HTML>
<p>
This is a simplified model of a tank. 
The top part is open to the environment at the fixed pressure 
<tt>p_ambient</tt>. Heat transfer to the environment and to 
the tank walls is neglected.
The tank is filled with a single or multiple-substance liquid, 
assumed to have uniform temperature and mass fractions.</p> 
<p>
Inlet and outlet connections are situated at the bottom of the tank. The following assumptions are made:
</p>
<ul>
<li>Incompressible medium, liquid fluid is of uniform density and uniform temperature</li>
<li>Heat transfer to the environment is neglected</li>
<li>Kinetic energy of the fluid in the tank is neglected, the cross sectional area of the tank is larger than the cross sectional area of the inlet/outlet by several orders of magnitude</li>
<li>No air is leaving the tank through the ports, if the liquid level drops below zero the simulation stops.</li>
<li>No liquid is leaving the tank through the open top, if the liquid level growth over the height the simulation stops.</li>
</ul>
<p>By default the port pressure is the pressure just after the outlet (or just before the inlet) in the attached pipe. The hydraulic resistances <tt>zeta_in</tt> and <tt>zeta_out</tt> determine the dissipative pressure drop between tank and port depending on the direction of mass flow. The default values (zeta_in=1, zeta_out=0) assume no dissipation at the tank outlet (ideal smooth opening) and total dissipation of kinetic energy at the tank inlet. Larger values are found for sharp edged openings and non-uniform velocity distributions in the pipe. A large selection of possible cases are listed in <i>[Idelchik, Handbook of Hydraulic Resistance, 2004]</i>. If the flag <tt>static_pressure_at_port</tt> in the <tt>Advanced</tt> menu is set to true, the port pressure represents the static head at the bottom of the tank. The relationship between pressure drop and mass flow rate must then be provided by the connected component. </p>  
 
 
</HTML>", revisions="<html>
<ul>
<li><i>Jan. 6, 2006</i> by Katja Poschlad, Manuel Remelhe (AST Uni Dortmund), 
   Martin Otter (DLR):<br> 
   Implementation based on former tank model.</li>
<li><i>Apr. 25, 2006</i> by Katrin Pr&ouml;l&szlig; (TUHH):<br>
Limitation to bottom ports only, added inlet and outlet loss factors.</li>
<li><i>Oct. 29, 2007</i> by Carsten Heinrich (ILK Dresden):<br>
Adapted to the new fluid library interfaces: 
<ul> <li>FluidPorts_b is used instead of FluidPort_b (due to it is defined as an array of ports)</li>
    <li>Port name changed from port to ports</li></ul>Updated documentation.</li>
 
</ul>
</html>"),
      Diagram,
      uses(Modelica(version="2.2.1"), Modelica_Fluid(version="0.952")),
      Coordsys(grid=[1,1], scale=0.2));
end OpenTank;
  
model Tank 
    "Open tank with top and bottom inlet/outlet ports at a defineable height" 
    
  import Modelica.Constants;
  import Modelica_Fluid.PressureLosses.BaseClasses.lossConstant_D_zeta;
  import Modelica_Fluid.Utilities.regRoot2;
  import Modelica_Fluid.Volumes.BaseClasses.TankPortData;
    
  replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
    annotation (choicesAllMatching=true);
    
  SI.Height level(stateSelect=StateSelect.prefer, start=level_start) 
      "Fluid level in the tank";
    
//Tank geometry  
    parameter SI.Height levelMax "Maximum level of tank before it overflows";
    parameter SI.Area area "Area of tank";
    parameter SI.Volume V0=0 "Volume of the liquid when level = 0";
    
//Port definitions 
    parameter Integer nTopPorts(min=1) = 1 
      "Number of inlet ports above levelMax (>= 1)";
    
    Modelica_Fluid.Interfaces.FluidPorts_a topPorts[nTopPorts](
    redeclare package Medium = Medium,
    m_flow(each start=0, each min=0),
    mXi_flow(each start=0, each min=0)) 
      "Inlet ports over levelMax at top of tank (fluid flows only from the port in to the tank)"
    annotation (extent=[-5,85; 5,125],  rotation=90);
    
    parameter Modelica_Fluid.Volumes.BaseClasses.TankPortData portsData[:] = {TankPortData(diameter=0)} 
      "Data of inlet/outlet ports at side and bottom of tank";
    
    Modelica_Fluid.Interfaces.FluidPorts_b ports[size(portsData,1)](
    redeclare package Medium = Medium,
    m_flow(each start=0),
    mXi_flow(each start=0)) 
      "inlet/outlet ports at bottom or side of tank (fluid flows in to or out of port; a port might be above the fluid level)"
    annotation (extent=[-5,-125; 5,-85],  rotation=90);
    
//Ambient  
   outer Modelica_Fluid.Ambient ambient "Ambient conditions";
   parameter Medium.AbsolutePressure p_ambient=ambient.default_p_ambient 
      "Tank surface pressure" 
    annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));
   parameter Medium.Temperature T_ambient=ambient.default_T_ambient 
      "Tank surface Temperature" 
    annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));
    
//Initialization
    parameter Types.Init.Temp initType=Types.Init.InitialValues 
      "Initialization option" 
    annotation(Evaluate=true,Dialog(tab = "Ambient and Initialization", group = "Initialization"));
    parameter SI.Height level_start(min=0) "Start value of tank level" 
    annotation(Dialog(tab="Ambient and Initialization", group = "Initialization"));
    parameter Medium.Temperature T_start=T_ambient "Start value of temperature"
    annotation(Dialog(tab = "Ambient and Initialization", group = "Initialization"));
    parameter Medium.MassFraction X_start[Medium.nX]=Medium.X_default 
      "Start value of mass fractions m_i/m" 
    annotation (Dialog(tab="Ambient and Initialization", group = "Initialization", enable=Medium.nXi > 0));
    
// Advanced  
    parameter Real hysteresisFactor(min=0) = 0.1 
      "Hysteresis for empty pipe = diameter*hysteresisFactor" 
    annotation(Dialog(tab="Advanced", group="Numerical properties"));
    parameter SI.MassFlowRate m_flow_small(min=0) = 1e-5 
      "Regularization range at zero mass flow rate" 
    annotation(Dialog(tab="Advanced", group="Numerical properties"));
    parameter Boolean stiffCharacteristicForEmptyPort = true 
      "=true, if steep pressure loss characteristic for empty pipe port" 
    annotation(Dialog(tab="Advanced", group="Numerical properties"), Evaluate=true);
    parameter Real zetaLarge(min=0) = 1e5 
      "Large pressure loss factor if mass flows out of empty pipe port" 
    annotation(Dialog(tab="Advanced", group="Numerical properties", enable=stiffCharacteristicForEmptyPort));
    
//Tank properties  
     final parameter Integer nPorts = size(ports,1) 
      "Number of inlet/outlet ports";
     final parameter Medium.SpecificEnthalpy h_start=Medium.specificEnthalpy_pTX(
        p_ambient,
        T_start,
        X_start) annotation(Hide=true);
    Medium.BaseProperties medium(
      preferredMediumStates=true,
      p(start=p_ambient),
      T(start=T_start),
      h(start=h_start),
      Xi(start=X_start[1:Medium.nXi]));
    SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";
    SI.Energy U(stateSelect=StateSelect.never) "Internal energy of tank volume";
    SI.Mass m(stateSelect=StateSelect.never) "Mass of fluid in tank";
    SI.Mass mXi[Medium.nXi](each stateSelect=StateSelect.never) 
      "Masses of independent components in the fluid";
    
  protected 
    parameter SI.Area bottomArea[nPorts]=Constants.pi*{(portsData[i].diameter/2)^2 for i in 1:nPorts};
    parameter SI.Diameter ports_emptyPipeHysteresis[nPorts] = portsData.diameter*hysteresisFactor;
    SI.Length levelAbovePort[nPorts] "Height of fluid over bottom ports";
    Boolean ports_m_flow_out[nPorts](each start = true, each fixed=true);
    Boolean aboveLevel[nPorts] "= true, if level >= ports[i].portLevel";
    Real zeta_out[nPorts];
equation 
  assert(level <= levelMax, "Tank starts to overflow (level = levelMax = " + String(level) + ")");
  assert(m>=0, "Mass in tank is zero");
    
  // Total quantities
    medium.p = p_ambient;
    V = area*level + V0 "Volume of fluid";
    m = V*medium.d "Mass of fluid";
    mXi = m*medium.Xi "Mass of fluid components";
    U = m*medium.u "Internal energy of fluid";
    
  // Mass balances
    der(m) = sum(topPorts.m_flow) + sum(ports.m_flow);
    for i in 1:Medium.nXi loop
      der(mXi[i]) = sum(topPorts.mXi_flow[i]) + sum(ports.mXi_flow[i]);
    end for;
    
  // Energy balance
    if Medium.singleState then
      der(U) = sum(topPorts.H_flow) + sum(ports.H_flow);
                               //Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)
    else
      der(U) = sum(topPorts.H_flow) + sum(ports.H_flow) - p_ambient*der(V);
    end if;
    
  // Properties at top ports
    for i in 1:nTopPorts loop
       // It is assumed that fluid flows only into one of the top ports and never out of it 
       topPorts[i].H_flow   = semiLinear(topPorts[i].m_flow, topPorts[i].h, h_start);
       topPorts[i].mXi_flow = semiLinear(topPorts[i].m_flow, topPorts[i].Xi, X_start[1:Medium.nXi]);
       topPorts[i].p        = p_ambient;
/*
       assert(topPorts[i].m_flow > -1, "Mass flows out of tank via topPorts[" + String(i) + "]\n" +
                                         "This indicates a wrong model");
*/
    end for;
    
  // Properties at bottom ports
    for i in 1:nPorts loop
       ports[i].H_flow = semiLinear(ports[i].m_flow, ports[i].h, medium.h);
       ports[i].mXi_flow = semiLinear(ports[i].m_flow, ports[i].Xi, medium.Xi);
       aboveLevel[i] = level >= (portsData[i].portLevel + ports_emptyPipeHysteresis[i])
                       or pre(aboveLevel[i]) and level >= (portsData[i].portLevel - ports_emptyPipeHysteresis[i]);
       levelAbovePort[i] = if aboveLevel[i] then level - portsData[i].portLevel else 0;
      
       if stiffCharacteristicForEmptyPort then
          // If port is above fluid level, use large zeta if fluid flows out of port (= small mass flow rate)
          zeta_out[i] = 1 + (if aboveLevel[i] then 0 else zetaLarge);
          ports[i].p = p_ambient + levelAbovePort[i]*ambient.g*medium.d
                               + Modelica_Fluid.Utilities.regSquare2(ports[i].m_flow, m_flow_small,
                                     lossConstant_D_zeta(portsData[i].diameter, 0.01)/medium.d,
                                     lossConstant_D_zeta(portsData[i].diameter, zeta_out[i])/medium.d);
          ports_m_flow_out[i] = false;
        
       else
          // Handling according to Remelhe/Poschlad
          ports_m_flow_out[i] = (pre(ports_m_flow_out[i]) and not ports[i].p>p_ambient)
                                     or ports[i].m_flow < -1e-6;
         if aboveLevel[i] then
             ports[i].p = p_ambient + levelAbovePort[i]*ambient.g*medium.d -
                               smooth(2,noEvent(if ports[i].m_flow < 0 then ports[i].m_flow^2/
                                     (2*medium.d*bottomArea[i]^2) else 0));
         else
            if pre(ports_m_flow_out[i]) then
               ports[i].m_flow = 0;
            else
               ports[i].p = p_ambient;
            end if;
         end if;
          zeta_out[i] =0;
       end if;
     end for;
    
initial equation 
    for i in 1:nPorts loop
       pre(aboveLevel[i]) = level_start >= portsData[i].portLevel;
    end for;
    
    if initType == Types.Init.NoInit then
    // no initial equations
    elseif initType == Types.Init.InitialValues then
      level = level_start;
      medium.T = T_start;
      medium.Xi = X_start[1:Medium.nXi];
    elseif initType == Types.Init.SteadyState then
      der(level) = 0;
      der(medium.T) = 0;
      der(medium.Xi) = zeros(Medium.nXi);
    elseif initType == Types.Init.SteadyStateHydraulic then
      der(level) = 0;
      medium.T = T_start;
      medium.Xi = X_start[1:Medium.nXi];
    else
      assert(false, "Unsupported initialization option");
    end if;
    
    annotation (
      Icon(
        Rectangle(extent=[-100,-100; 100,100], style(
            color=7,
            rgbcolor={255,255,255},
            fillColor=7,
            rgbfillColor={255,255,255})),
          Rectangle(extent=DynamicSelect([-100,-100; 100,0], [-100,-100; 100,(-100
                 + 200*level/levelMax)]), style(
              color=69,
              rgbcolor={0,127,255},
              fillColor=71,
              rgbfillColor={85,170,255},
              fillPattern=1)),
        Text(
          extent=[-94,19; 96,-1],
        string=DynamicSelect(" ", realString(level, 1, 3)),
        style(color=0, rgbcolor={0,0,0})),
        Line(points=[-100,100; 100,100], style(
            color=0,
            rgbcolor={0,0,0},
            pattern=3)),
      Text(
          extent=[-94,90;95,60],
          style(color=3, rgbcolor={0,0,255}),
          string="%name"),
      Text(
        extent=[-95,-85; 95,-65],
        style(color=0),
          string="%level_start"),
      Text(
        extent=[-95,-55; 95,-35],
        style(color=0),
          string="level_start ="),
      Text(extent=[-95,50; 95,30], string="level =",
        style(color=0, rgbcolor={0,0,0})),
        Line(points=[-100,100; -100,-100; 100,-100; 100,100], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=69,
            rgbfillColor={0,127,255},
            fillPattern=1))),
      Documentation(info="<HTML>
<p> 
Model of a tank that is open to the environment at the fixed pressure
<tt>p_ambient</tt>. Heat transfer to the environment and to 
the tank walls is neglected.
The tank is filled with a single or multiple-substance liquid, 
assumed to have uniform temperature and mass fractions.
</p> 
 
<p>
At the top of the tank over the maximal fill level <b>levelMax</b> 
a vector of FluidPorts, called <b>topPorts</b>, is present.
The assumption is made that fluid flows always in to the tank via these
ports (and never back in to the connector).
If the tank has no top ports, set <b>nTopPorts</b> = 1, and do not
connect to this port (the default connection semantics of Modelica
leads to a behaviour as if the port would not be present; 
the reason is that some tools do currently no support zero sized
connectors).
</p>
 
<p>
The vector of connectors <b>ports</b> are fluid ports at the bottom
and side of the tank at a defineable height. Fluid can flow either out
of or in to this port. The fluid level of the tank may be below
one of these ports. This case is approximated by introducing a
large pressure flow coefficient so that the mass flow rate
through this port is very small in this case.
</p>
 
<p>
If the tank starts to over flow (i.e., level > levelMax), an
assertion is triggered.
</p>
 
<p>
When the diagram layer is open in the plot environment, the
level of the tank is dynamically visualized. Note, the speed
of the diagram animation in Dymola can be set via command
<b>animationSpeed</b>(), e.g., animationSpeed(speed = 10)
</p>
</HTML>", revisions="<html>
<ul>
<li><i>Jul. 29, 2006</i> by Martin Otter (DLR):<br> 
   Improved handling of ports that are above the fluid level and
   simpler implementation.</li>
 
<li><i>Jan. 6, 2006</i> by Katja Poschlad, Manuel Remelhe (AST Uni Dortmund), 
   Martin Otter (DLR):<br> 
   Implementation based on former tank model but with several improvements
   (top, bottom, side ports; correctly treating kinetic energy for outlet
   and total dissipation for inlet; ports can be above the fluid level).</li>
</ul>
</html>"),
      Diagram,
      uses(Modelica(version="2.2.1"), Modelica_Fluid(version="0.952")),
      Coordsys(grid=[1,1], scale=0.2));
end Tank;
  
  package BaseClasses 
    extends Modelica_Fluid.Icons.BaseClassLibrary;
    
    record TankPortData 
      "Data to describe inlet/outlet pipes at the bottom or side of the tank" 
          extends Modelica.Icons.Record;
      
      parameter SI.Diameter diameter 
        "Inner (hydraulic) diameter of inlet/outlet port";
      parameter SI.Height portLevel=0 
        "level of inlet/outlet port (height over the tank base)";
    end TankPortData;
    
      partial model PartialLumpedVolume 
      "Mixing volume with inlet and outlet ports (flow reversal is allowed)" 
      import Modelica_Fluid.Types;
        replaceable package Medium = 
          Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
            annotation (choicesAllMatching = true);
        parameter Types.Init.Temp initType=
                  Types.Init.NoInit "Initialization option" 
          annotation(Evaluate=true, Dialog(tab = "Initialization"));
        parameter Medium.AbsolutePressure p_start = Medium.p_default 
        "Start value of pressure" 
          annotation(Dialog(tab = "Initialization"));
        parameter Boolean use_T_start = true 
        "= true, use T_start, otherwise h_start" 
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
      
        parameter Types.FlowDirection.Temp flowDirection=Types.FlowDirection.
            Bidirectional 
        "Unidirectional (port_a -> port_b) or bidirectional flow component" 
         annotation(Dialog(tab="Advanced"));
         parameter Boolean allowFlowReversal=flowDirection == Types.FlowDirection.
            Bidirectional 
        "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
         annotation(Evaluate=true, Hide=true);
      
        Interfaces.FluidStatePort_a port_a(
                                      redeclare package Medium = Medium, m_flow(min=
                if allowFlowReversal then -Modelica.Constants.inf else 0)) 
        "Fluid inlet port" annotation (extent=[-112,-10; -92,10]);
        Interfaces.FluidStatePort_b port_b(
                                      redeclare package Medium = Medium, m_flow(max=
                if allowFlowReversal then +Modelica.Constants.inf else 0)) 
        "Fluid outlet port" annotation (extent=[90,-10; 110,10]);
        Medium.BaseProperties medium(
          preferredMediumStates=true,
          p(start=p_start),
          h(start=h_start),
          T(start=T_start),
          Xi(start=X_start[1:Medium.nXi]));
        SI.Energy U "Internal energy of fluid";
        SI.Mass m "Mass of fluid";
        SI.Mass mXi[Medium.nXi] "Masses of independent components in the fluid";
        SI.Volume V_lumped "Volume";
      
    protected 
        SI.HeatFlowRate Qs_flow 
        "Heat flow across boundaries or energy source/sink";
        SI.Power Ws_flow "Work flow across boundaries or source term";
        annotation (
          Icon(Text(extent=[-144,178; 146,116], string="%name"), Text(
              extent=[-130,-108; 144,-150],
              style(color=0),
              string="V=%V")),
          Documentation(info="<html>
Base class for an ideally mixed fluid volume with two ports and the ability to store mass and energy. The following source terms are part of the energy balance and must be specified in the extending class:
<ul>
<li><tt>Qs_flow</tt>, e.g. convective or latent heat flow rate across segment boundary, and</li> <li><tt>Ws_flow</tt>, work term, e.g. p*der(V) if the volume is not constant</li>
</ul>
The component volume <tt>V_lumped</tt> is also a variable which needs to be set in the extending class to complete the model.
</html>"),Diagram);
      
      equation 
      // boundary conditions
        port_a.p = medium.p;
        port_b.p = medium.p;
        port_a.H_flow = semiLinear(
          port_a.m_flow,
          port_a.h,
          medium.h);
        port_b.H_flow = semiLinear(
          port_b.m_flow,
          port_b.h,
          medium.h);
        port_a.mXi_flow = semiLinear(
          port_a.m_flow,
          port_a.Xi,
          medium.Xi);
        port_a.mC_flow = semiLinear(
          port_a.m_flow,
          port_a.C,
          port_b.C);
        port_b.mXi_flow = semiLinear(
          port_b.m_flow,
          port_b.Xi,
          medium.Xi);
      
      // Total quantities
        m = V_lumped*medium.d;
        mXi = m*medium.Xi;
        U = m*medium.u;
      
      // Mass and energy balance
        der(m) = port_a.m_flow + port_b.m_flow;
        der(mXi) = port_a.mXi_flow + port_b.mXi_flow;
        der(U) = port_a.H_flow + port_b.H_flow + Qs_flow + Ws_flow;
        zeros(Medium.nC) = port_a.mC_flow+port_b.mC_flow;
      
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
      end PartialLumpedVolume;
  end BaseClasses;
  annotation (Documentation(info="<html>
 
</html>"));
end Volumes;
