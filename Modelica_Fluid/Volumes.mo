package Volumes 
   extends Modelica_Fluid.Icons.VariantLibrary;
  
    model MixingVolume 
    "Mixing volume with inlet and outlet ports (flow reversal is allowed)" 
    extends Modelica_Fluid.Interfaces.PartialLumpedVolume(
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
</html>"),
      Diagram);
    equation 
    thermalPort.T = medium.T;
    Qs_flow = thermalPort.Q_flow;
    end MixingVolume;
  
model OpenTank "Open tank with inlet/outlet ports at the bottom" 
    import SI = Modelica.SIunits;
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
    Modelica_Fluid.Interfaces.FluidPort_a port[n_ports](
      redeclare package Medium = Medium,
      m_flow(each start=0),
      mXi_flow(each start=0)) 
    annotation (extent=[-12,-109; 8,-89]);
    
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
    der(m) = sum(port.m_flow);
    for i in 1:Medium.nXi loop
      der(mXi[i]) = sum(port[:].mXi_flow[i]);
    end for;
    
  // Energy balance
    if Medium.singleState then
      der(U) = sum(port.H_flow);
                               //Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)
    else
      der(U) = sum(port.H_flow) - p_ambient*der(V);
    end if;
    assert(level <= height, "Tank is full (level = height = " + String(level) + ")");
    assert(level > 0, "Tank is empty (level = 0), tank model is not designed to allow air flow through ports");
    
//Determine port properties  
    p_static = level*ambient.g*medium.d + p_ambient;
    for i in 1:n_ports loop
      port[i].H_flow = semiLinear(
        port[i].m_flow,
        port[i].h,
        medium.h);
      port[i].mXi_flow = semiLinear(
        port[i].m_flow,
        port[i].Xi,
        medium.Xi);
      if p_static_at_port then
        port[i].p = p_static;
      else
       port[i].p = p_static - smooth(2, noEvent(port[i].m_flow^2/(2*medium.d*
          pipeArea[i]^2)*(if port[i].m_flow < 0 then (1 + zeta_out[i]) else (1
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
      defaultComponentName="tank",
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
</ul>
<p>By default the port pressure is the pressure just after the outlet (or just before the inlet) in the attached pipe. The hydraulic resistances <tt>zeta_in</tt> and <tt>zeta_out</tt> determine the dissipative pressure drop between tank and port depending on the direction of mass flow. The default values (zeta_in=1, zeta_out=0) assume no dissipation at the tank outlet (ideal smooth opening) and total dissipation of kinetic energy at the tank inlet. Larger values are found for sharp edged openings and non-uniform velocity distributions in the pipe. A large selection of possible cases are listed in <i>[Idelchik, Handbook of Hydraulic Resistance, 2004]</i>. If the flag <tt>static_pressure_at_port</tt> in the <tt>Advanced</tt> menu is set to true, the port pressure represents the static head at the bottom of the tank. The relationship between pressure drop and mass flow rate must then be provided by the connected component. </p>  
 
 
</HTML>", revisions="<html>
<ul>
<li><i>Jan. 6, 2006</i> by Katja Poschlad, Manuel Remelhe (AST Uni Dortmund), 
   Martin Otter (DLR):<br> 
   Implementation based on former tank model.</li>
<li><i>Apr. 25, 2006</i> by Katrin Pr&ouml;l&szlig; (TUHH):<br>
Limitation to bottom ports only, added inlet and outlet loss factors.</li>
</ul>
</html>"),
      Diagram,
      uses(Modelica(version="2.2.1"), Modelica_Fluid(version="0.952")),
      Coordsys(grid=[1,1], scale=0.2));
end OpenTank;
  
model OpenTank2 
    "Open tank with top and bottom inlet/outlet ports at a defineable height" 
  import SI = Modelica.SIunits;
  import Modelica.Constants;
  replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
    annotation (choicesAllMatching=true);
    
//Tank geometry  
    parameter SI.Area area "Area of tank";
    parameter SI.Volume V0=0 "Volume of the liquid when level = 0";
    parameter SI.Height levelMax "Maximum level of tank before it overflows";
    
//Port definitions 
  /*
     {Modelica_Fluid.Volumes.BaseClasses.TankTopPortData()} 
 
           {Modelica_Fluid.Volumes.BaseClasses.TankBottomPortData()} 
 */
    parameter Modelica_Fluid.Volumes.BaseClasses.TankTopPortData topPortData[:] 
      "Data of inlet ports at top of tank";
    parameter Modelica_Fluid.Volumes.BaseClasses.TankBottomPortData 
      bottomPortData[
                   :] "Data of inlet/outlet ports at side and bottom of tank";
    
    Modelica_Fluid.Interfaces.FluidPort_a topPort[size(topPortData,1)](
      redeclare package Medium = Medium,
      m_flow(each start=0),
      mXi_flow(each start=0)) 
    annotation (extent=[-10,90; 10,110]);
    
    Modelica_Fluid.Interfaces.FluidPort_b bottomPort[size(bottomPortData,1)](
      redeclare package Medium = Medium,
      m_flow(each start=0),
      mXi_flow(each start=0)) 
    annotation (extent=[-10,-110; 10,-90]);
    
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
    parameter SI.Height level_start "Start value of tank level" 
    annotation(Dialog(tab="Ambient and Initialization", group = "Initialization"));
    parameter Medium.Temperature T_start=T_ambient "Start value of temperature"
    annotation(Dialog(tab = "Ambient and Initialization", group = "Initialization"));
    parameter Medium.MassFraction X_start[Medium.nX]=Medium.X_default 
      "Start value of mass fractions m_i/m" 
    annotation (Dialog(tab="Ambient and Initialization", group = "Initialization", enable=Medium.nXi > 0));
    
//Transformation of kinetic energy
    parameter Boolean p_static_at_port=false 
      "= true if kinetic energy and dissipation is neglected at inlet/outlet ports"
                                                                                   annotation(Evaluate=true, Dialog(tab="Advanced"));
    
//Tank properties  
    final parameter Integer nTop = size(topPortData,1);
    final parameter Integer nBottom = size(bottomPortData,1);
    
    final parameter Medium.SpecificEnthalpy h_start=Medium.specificEnthalpy_pTX(
        p_ambient,
        T_start,
        X_start);
    Medium.BaseProperties medium(
      preferredMediumStates=true,
      p(start=p_ambient),
      T(start=T_start),
      h(start=h_start),
      Xi(start=X_start[1:Medium.nXi]));
    SI.Height level(stateSelect=StateSelect.prefer, start=level_start) 
      "Level height of tank";
    SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";
    SI.Energy U(stateSelect=StateSelect.never) "Internal energy of tank volume";
    SI.Mass m(stateSelect=StateSelect.never) "Mass of fluid in tank";
    SI.Mass mXi[Medium.nXi](each stateSelect=StateSelect.never) 
      "Masses of independent components in the fluid";
    SI.Length levelAbovePort[nBottom] "Height of fluid over bottom ports";
    Medium.MassFlowRate overflow_m_flow 
      "Mass flow rate that overflows thank if level >= levelMax";
    Medium.MassFlowRate overflow_mXi_flow[Medium.nXi] 
      "Substance mass flow rate that overflows";
    Medium.EnthalpyFlowRate overflow_H_flow 
      "Enthalpy flow rate of overflowing fluid";
    Boolean bottomPort_m_flow_out[nBottom];
    Medium.AbsolutePressure p_static[nBottom];
    Boolean aboveLevel[nBottom] "= true, if level >= bottomPort[i].level";
    Real bottomPort_s[nBottom];
  protected 
    parameter SI.Area topArea[nTop]=Constants.pi*{(topPortData[i].diameter/2)^2 for i in 1:nTop};
    parameter SI.Area bottomArea[nTop]=Constants.pi*{(bottomPortData[i].diameter/2)^2 for i in 1:nBottom};
    parameter Real k_small(min=0) = 1e-5 
      "Small regularization range if tank level is below bottom_height or side_height; k_small = 0 gives ideal switch";
    
equation 
  //Total quantities
    medium.p = p_ambient;
    V = area*level + V0 "Volume of fluid";
    m = V*medium.d "Mass of fluid";
    mXi = m*medium.Xi "Mass of fluid components";
    U = m*medium.u "Internal energy of fluid";
    
    overflow_m_flow = 0;
    overflow_mXi_flow = zeros(Medium.nXi);
    overflow_H_flow = 0;
    
  // Mass balances
    der(m) = sum(topPort.m_flow) + sum(bottomPort.m_flow) + overflow_m_flow;
    for i in 1:Medium.nXi loop
      der(mXi[i]) = sum(topPort.mXi_flow[i]) + sum(bottomPort.mXi_flow[i]) + overflow_mXi_flow[i];
    end for;
    
  // Energy balance
    if Medium.singleState then
      der(U) = sum(topPort.H_flow) + sum(bottomPort.H_flow) + overflow_H_flow;
                               //Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)
    else
      der(U) = sum(topPort.H_flow) + sum(bottomPort.H_flow) + overflow_H_flow - p_ambient*der(V);
    end if;
    
  // Properties at top ports
    for i in 1:nTop loop
       /* It is assumed that fluid flows only into one of the top ports and never out of it 
       */
       topPort[i].H_flow   = semiLinear(topPort[i].m_flow, topPort[i].h, h_start);
       topPort[i].mXi_flow = semiLinear(topPort[i].m_flow, topPort[i].Xi, X_start[1:Medium.nXi]);
       topPort[i].p        = p_ambient;
       assert(topPort[i].m_flow > -1e-3, "Mass flows out of tank via topPort[" + String(i) + "]\n" +
                                         "This indicates a wrong model");
    end for;
    
  // Properties at bottom ports  
    for i in 1:nBottom loop
       bottomPort[i].H_flow = semiLinear(bottomPort[i].m_flow, bottomPort[i].h, medium.h);
       bottomPort[i].mXi_flow = semiLinear(bottomPort[i].m_flow, bottomPort[i].Xi, medium.Xi);
       levelAbovePort[i] = level - bottomPortData[i].level;
       aboveLevel[i] = levelAbovePort[i] >= 0;
      
       if p_static_at_port then
          assert(aboveLevel[i], "Level (= " + String(level) + ") is below bottomPort[" + String(i) + "].level.\n" +
                                "This is not allowed if p_static_at_port = true");
          p_static[i] = p_ambient + levelAbovePort[i]*ambient.g*medium.d;
          bottomPort_m_flow_out[i] = bottomPort[i].m_flow < 0;
          bottomPort[i].p = p_static[i];
          bottomPort_s[i] = 0;
       else
          p_static[i] = p_ambient + (if aboveLevel[i] then levelAbovePort[i]*ambient.g*medium.d else 0);
          bottomPort_m_flow_out[i] = bottomPort[i].m_flow < 0;
          bottomPort[i].p = 0;
/*
          bottomPort[i].p = p_static - smooth(2, noEvent(bottomPort[i].m_flow^2/(2*medium.d*bottomArea[i]^2)*
                                       (if  bottomPort_m_flow_out[i] then bottomPortData[i].zeta_out else 
                                                                          bottomPortData[i].zeta_in)));
*/
          bottomPort_s[i] = 0;
        
/*
          bottomPort_m_flow_out[i] = bottomPort_s[i] <= 0;
          bottomPort[i].m_flow = if aboveLevel[i] then bottomPort_s[i] else 
                                    (if bottomPort_m_flow_out[i] then k_small else 1)*bottomPort_s[i];
          bottomPort[i].p = p_static[i] +
                           (if aboveLevel[i] then 
                               -smooth(2,if bottomPort_m_flow_out[i] then 
                                             bottomPort_s[i]*abs(bottomPort_s[i])/(2*medium.d*bottomArea[i]^2) else 
                                             k_small*bottomPort_s[i]) else 
                                 (if bottomPort_m_flow_out[i] then bottomPort_s[i] else k_small*bottomPort_s[i]));
*/
       end if;
    end for;
    
/*
   OpenTank (Katrin):
       if p_static_at_port then
          bottomPort[i].p = p_static[i];
       else
          port[i].p = p_static - smooth(2, noEvent(port[i].m_flow^2/(2*medium.d*
             pipeArea[i]^2)*(if port[i].m_flow < 0 then (1 + zeta_out[i]) else (1
             - zeta_in[i]))));
       end if;
 
   TankAttachment (Martin):
       m_flow_out = s <= 0;
       if aboveLevel >= 0 then
          m_flow = m_flow_nominal*s "equation to compute s, which is a dummy in this branch";
          port.p - p_ambient = aboveLevel*fluidOptions.g*d  -
                               smooth(2,if m_flow_out then s*abs(s)*m_flow_nominal^2/(2*d*pipeArea^2) else k_small*m_flow_nominal*s);
       else
          m_flow = (if m_flow_out then k_small*p_nominal else m_flow_nominal)*s;
          port.p - p_ambient = (if m_flow_out then p_nominal else k_small*m_flow_nominal)*s;
       end if;
*/
    
initial equation 
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
      defaultComponentName="tank",
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
          extent=[-95,84; 95,54],
          string="%name",
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=69,
            rgbfillColor={0,127,255},
            fillPattern=1)),
        Text(
          extent=[-129,47; 130,33],
          style(color=0),
          Line(points=[-100,100; 100,100], style(
              color=0,
              rgbcolor={0,0,0},
              pattern=3)),
          string="start = %level_start m"),
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
</ul>
<p>By default the port pressure is the pressure just after the outlet (or just before the inlet) in the attached pipe. The hydraulic resistances <tt>zeta_in</tt> and <tt>zeta_out</tt> determine the dissipative pressure drop between tank and port depending on the direction of mass flow. The default values (zeta_in=1, zeta_out=0) assume no dissipation at the tank outlet (ideal smooth opening) and total dissipation of kinetic energy at the tank inlet. Larger values are found for sharp edged openings and non-uniform velocity distributions in the pipe. A large selection of possible cases are listed in <i>[Idelchik, Handbook of Hydraulic Resistance, 2004]</i>. If the flag <tt>static_pressure_at_port</tt> in the <tt>Advanced</tt> menu is set to true, the port pressure represents the static head at the bottom of the tank. The relationship between pressure drop and mass flow rate must then be provided by the connected component. </p>  
 
 
</HTML>", revisions="<html>
<ul>
<li><i>Jan. 6, 2006</i> by Katja Poschlad, Manuel Remelhe (AST Uni Dortmund), 
   Martin Otter (DLR):<br> 
   Implementation based on former tank model.</li>
<li><i>Apr. 25, 2006</i> by Katrin Pr&ouml;l&szlig; (TUHH):<br>
Limitation to bottom ports only, added inlet and outlet loss factors.</li>
</ul>
</html>"),
      Diagram,
      uses(Modelica(version="2.2.1"), Modelica_Fluid(version="0.952")),
      Coordsys(grid=[1,1], scale=0.2));
end OpenTank2;
  
  package BaseClasses 
    extends Modelica_Fluid.Icons.BaseClassLibrary;
    record TankTopPortData "Data to describe inlet pipe at top of tank" 
      import SI = Modelica.SIunits;
      extends Modelica.Icons.Record;
      
      parameter SI.Diameter diameter "Inner (hydraulic) diameter of inlet port";
      parameter SI.Height level 
        "Level of inlet port (height over the tank base; not yet used)";
    end TankTopPortData;
    
    record TankBottomPortData 
      "Data to describe inlet/outlet pipe of tank side or tank bottom" 
      import SI = Modelica.SIunits;
      extends Modelica.Icons.Record;
      
      parameter SI.Diameter diameter 
        "Inner (hydraulic) diameter of inlet/outlet port";
      parameter SI.Height level 
        "level of inlet/outlet port (height over the tank base)";
      parameter Real zeta_in=1 
        "<html>Hydraulic pressure loss factor zeta if fluid flows <b>in to</b> tank<br>(= 1 for total dissipation of kinetic energy and uniform flow distribution in pipe)</html>";
      parameter Real zeta_out=0 
        "<html>Hydraulic pressure loss factor zeta if fluid flows <b>out of</b> tank<br>(= 0 for ideal smooth outlet)</html>";
    end TankBottomPortData;
  end BaseClasses;
end Volumes;
