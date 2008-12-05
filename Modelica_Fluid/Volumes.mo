within Modelica_Fluid;
package Volumes "Generic volume, tank and other volume type components"
   extends Modelica_Fluid.Icons.VariantLibrary;

    model Volume "Fixed volume with ports"
      extends Modelica_Fluid.Volumes.BaseClasses.PartialLumpedVolumePorts(Qs_flow = heatPort.Q_flow);

      parameter SI.Volume V "Volume";

    //Port definitions
      parameter Integer nPorts(min=1)=1 "Number of ports" annotation(Dialog(group="Ports"));
      parameter SI.Diameter portDiameters[nPorts] = ones(nPorts)
      "Inner (hydraulic) diameters of ports (array)"   annotation(Dialog(group="Ports",enable= not neglectPortDiameters));

    //Transformation of kinetic energy
        parameter Boolean neglectPortDiameters=true
      "=true, kinetic energy and dissipation is accounted for in port pressure"
                                                                                                          annotation(Evaluate=true, Dialog(tab="Assumptions"));
        parameter Real[nPorts] zeta_in=fill(0, nPorts)
      "Hydraulic resistance into volume, 1 for total dissipation of kinetic energy and uniform flow distribution in pipe"
                                                                                                          annotation(Dialog(tab="Assumptions",enable= not neglectPortDiameters));
        parameter Real[nPorts] zeta_out=fill(1, nPorts)
      "Hydraulic resistance out of volume, 0 for ideal smooth outlet"   annotation(Dialog(tab="Assumptions",enable= not neglectPortDiameters));

      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort
      "Thermal port" 
        annotation (Placement(transformation(extent={{-20,88},{20,108}}, rotation=0)));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{
              100,100}}), graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillPattern=FillPattern.Sphere,
            fillColor={170,213,255})}),
      Documentation(info="<html>
Ideally mixed volume of constant size with two fluid ports and one medium model. The flow properties are computed from the upstream quantities, pressures are equal in both nodes and the medium model. Heat transfer through a thermal port is possible, it equals zero if the port remains unconnected. The thermal port temperature is equal to the medium temperature.
</html>"),
      Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{
              100,100}}),
              graphics));

    equation
      fluidVolume = V;
      Ws_flow = 0;
      heatPort.T = medium.T;
    end Volume;

  model SweptVolume
    "varying cylindric volume depending on the postition of the piston"
    extends BaseClasses.PartialLumpedVolumePorts(Qs_flow = 0);
    parameter SI.Area pistonCrossArea "cross sectional area of pistion";
    parameter SI.Volume clearance "remaining volume at zero piston stroke";

    annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
              -100},{100,100}}),
                        graphics),
                         Icon(coordinateSystem(preserveAspectRatio=true,
            extent={{-100,-100},{100,100}}), graphics={
          Rectangle(
            extent={{-44,36},{44,-90}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            lineThickness=1,
            fillColor={170,213,255},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-44,60},{-40,60},{-40,-32},{-44,-32},{-44,60}},
            lineColor={95,95,95},
            smooth=Smooth.None,
            fillColor={135,135,135},
            fillPattern=FillPattern.Backward),
          Polygon(
            points={{40,60},{44,60},{44,-34},{40,-34},{40,60}},
            lineColor={95,95,95},
            smooth=Smooth.None,
            fillColor={135,135,135},
            fillPattern=FillPattern.Backward),
          Rectangle(
            extent={{-40,40},{40,30}},
            lineColor={95,95,95},
            fillColor={135,135,135},
            fillPattern=FillPattern.Forward),
          Rectangle(
            extent={{-6,92},{6,40}},
            lineColor={95,95,95},
            fillColor={135,135,135},
            fillPattern=FillPattern.Forward),
          Polygon(
            points={{-40,-86},{40,-86},{40,70},{44,70},{44,-90},{-44,-90},{-44,
                70},{-40,70},{-40,-86}},
            lineColor={95,95,95},
            smooth=Smooth.None,
            fillColor={135,135,135},
            fillPattern=FillPattern.Backward)}),
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
      "translation flange for piston" annotation (Placement(transformation(
            extent={{-10,90},{10,110}},   rotation=0)));

  equation
    assert(flange.s >= 0, "Piston stroke (given by flange.s) must not be smaller than zero!");

    // volume size
    fluidVolume = clearance + flange.s * pistonCrossArea;

    flange.f = (medium.p - system.p_ambient) * pistonCrossArea;

    // energy balances
    Ws_flow = medium.p * pistonCrossArea * (-der(flange.s));
  end SweptVolume;

model OpenTank "Open tank with inlet/outlet ports at the bottom"
    replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
    annotation (choicesAllMatching=true);

  parameter Modelica_Fluid.Types.Dynamics dynamicsType=system.dynamicsType
      "Dynamics option" 
    annotation(Evaluate=true, Dialog(tab = "Assumptions"));

//Transformation of kinetic energy
    parameter Boolean neglectPortDiameters=false
      "=true, kinetic energy and dissipation is accounted for in port pressure"
                                                                                                        annotation(Evaluate=true, Dialog(tab="Assumptions"));
    parameter Real[nPorts] zeta_in=fill(0, nPorts)
      "Hydraulic resistance into tank, 1 for total dissipation of kinetic energy and uniform flow distribution in pipe"
                                                                                                     annotation(Dialog(tab="Assumptions",enable=neglectPortDiameters==false));
    parameter Real[nPorts] zeta_out=fill(1, nPorts)
      "Hydraulic resistance out of tank, 0 for ideal smooth outlet" 
                                                                  annotation(Dialog(tab="Assumptions",enable=neglectPortDiameters==false));

//Tank geometry
    parameter SI.Height height "Height of tank";
    parameter SI.Area crossArea "Area of tank";
    parameter SI.Volume V0=0 "Volume of the liquid when the level is zero";

//Port definitions
    parameter Integer nPorts(min=1) = 1 "Number of bottom ports (min=1)" 
     annotation(Dialog(group="bottomPorts (= pipes at bottom of tank; in and out flow of tank)"));
    parameter SI.Diameter portDiameters[nPorts]
      "Inner (hydraulic) diameters of bottom ports (array)" 
     annotation(Dialog(group="bottomPorts (= pipes at bottom of tank; in and out flow of tank)", enable=n_bottomPorts > 0));
    Interfaces.FluidPorts_b ports[nPorts](
    redeclare package Medium = Medium,
    m_flow(each start=0)) 
    annotation (Placement(transformation(
        extent={{0,-20},{-10,20}},
        rotation=90,
        origin={0,-100})));

//Ambient
    parameter Medium.AbsolutePressure p_ambient=system.p_ambient
      "Tank surface pressure" 
    annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));
    parameter Medium.Temperature T_ambient=system.T_ambient
      "Tank surface Temperature" 
    annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));
    outer Modelica_Fluid.System system "System properties";

//Initialization
    parameter Types.Init initType=Types.Init.InitialValues
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
    Medium.EnthalpyFlowRate H_flow[nPorts]
      "Enthalpy flow rates from the bottom ports in to the tank";
    Medium.MassFlowRate mXi_flow[nPorts, Medium.nXi]
      "Substance mass flow rates from the bottom ports in to the tank";

  protected
    parameter Medium.SpecificEnthalpy h_start=Medium.specificEnthalpy_pTX(
        p_ambient,
        T_start,
        X_start);
    parameter SI.Area[nPorts] portArea=Modelica.Constants.pi/4*{portDiameters[
        i]^2 for i in 1:nPorts};

equation
  // Only one connection allowed to a port to avoid unwanted ideal mixing
/*
for i in 1:nPorts loop
  assert(cardinality(ports[i]) <= 1,"
ports[" + String(i) + "] of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller.
");
end for;
*/

  //Total quantities
    V = crossArea*level + V0 "Volume of fluid";
    m = V*medium.d "Mass of fluid";
    mXi = m*medium.Xi "Mass of fluid components";
    U = m*medium.u "Internal energy of fluid";
    medium.p = p_ambient;

  // Mass balances
  if dynamicsType < Types.Dynamics.SteadyStateMass then
    der(m) = sum(ports.m_flow);
    for i in 1:Medium.nXi loop
      der(mXi[i]) = sum(mXi_flow[:,i]);
    end for;
  else
    0 = sum(ports.m_flow);
    for i in 1:Medium.nXi loop
      0 = sum(mXi_flow[:,i]);
    end for;
  end if;

  // Energy balance
  if dynamicsType < Types.Dynamics.SteadyState then
    if Medium.singleState then
      der(U) = sum(H_flow)
          "Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)";
    else
      der(U) = sum(H_flow) - p_ambient*der(V);
    end if;
  else
    if Medium.singleState then
      0 = sum(H_flow)
          "Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)";
    else
      0 = sum(H_flow) - p_ambient*der(V);
    end if;
  end if;
  assert(level <= height, "Tank is full (level = height = " + String(level) + ")");
  assert(level > 0, "Tank is empty (level = 0), tank model is not designed to allow air flow through ports");

//Determine port properties
    p_static = level*system.g*medium.d + p_ambient;
    for i in 1:nPorts loop
       H_flow[i]     = ports[i].m_flow*actualStream(ports[i].h_outflow);
       mXi_flow[i,:] = ports[i].m_flow*actualStream(ports[i].Xi_outflow);
       if neglectPortDiameters then
         ports[i].p = p_static;
       else
         ports[i].p = p_static - smooth(2, noEvent(ports[i].m_flow^2/(2*medium.d*
            portArea[i]^2)*(if ports[i].m_flow < 0 then (1 + zeta_out[i]) else (1
             - zeta_in[i]))));
       end if;
       ports[i].h_outflow = medium.h;
       ports[i].Xi_outflow = medium.Xi;
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
      Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1},
          initialScale=0.2), graphics={
          Rectangle(
            extent={{-100,100},{100,0}},
            lineColor={255,255,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent=DynamicSelect({{-100,-100},{100,10}}, {{-100,-100},{100,(-100
                 + 200*level/height)}}),
            lineColor={0,127,255},
            fillColor={85,170,255},
            fillPattern=FillPattern.Solid),
          Line(points={{-100,100},{-100,-100},{100,-100},{100,100}}, color={0,0,
                0}),
          Text(
            extent={{-95,90},{95,60}},
            lineColor={0,0,255},
            textString="%name"),
          Text(
            extent={{-129,53},{130,39}},
            lineColor={0,0,0},
            textString="start = %level_start m"),
          Text(
            extent={{-95,30},{95,5}},
            lineColor={0,0,0},
            textString=DynamicSelect(" ", realString(
                level,
                1,
                integer(precision)))),
          Line(
            points={{-100,100},{100,100}},
            color={0,0,0},
            pattern=LinePattern.Dot)}),
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
      Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1},
          initialScale=0.2), graphics),
      uses(Modelica(version="2.2.1"), Modelica_Fluid(version="0.952")));

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

  parameter Modelica_Fluid.Types.Dynamics dynamicsType=system.dynamicsType
      "Dynamics option" 
    annotation(Evaluate=true, Dialog(tab = "Assumptions"));

  SI.Height level(stateSelect=StateSelect.prefer, start=level_start)
      "Fluid level in the tank";

//Tank geometry
    parameter SI.Height levelMax "Maximum level of tank before it overflows";
    parameter SI.Area crossArea "Area of tank";
    parameter SI.Volume V0=0 "Volume of the liquid when level = 0";

//Port definitions
    parameter Integer nTopPorts(min=1) = 1
      "Number of inlet ports above levelMax (>= 1)";

    Modelica_Fluid.Interfaces.FluidPorts_a topPorts[nTopPorts](
    redeclare package Medium = Medium,
    m_flow(each start=0, each min=0))
      "Inlet ports over levelMax at top of tank (fluid flows only from the port in to the tank)"
    annotation (Placement(transformation(
        extent={{0,-20},{10,20}},
        rotation=90,
        origin={0,100})));

    parameter Modelica_Fluid.Volumes.BaseClasses.TankPortData portsData[:] = {TankPortData(diameter=0.0001)}
      "Data of inlet/outlet ports at side and bottom of tank";

    Modelica_Fluid.Interfaces.FluidPorts_b ports[size(portsData,1)](
    redeclare package Medium = Medium,
    m_flow(each start=0))
      "inlet/outlet ports at bottom or side of tank (fluid flows in to or out of port; a port might be above the fluid level)"
    annotation (Placement(transformation(
        extent={{0,-20},{-10,20}},
        rotation=90,
        origin={0,-100})));

//Ambient
   outer Modelica_Fluid.System system "System properties";
   parameter Medium.AbsolutePressure p_ambient=system.p_ambient
      "Tank surface pressure" 
    annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));
   parameter Medium.Temperature T_ambient=system.T_ambient
      "Tank surface Temperature" 
    annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));

//Initialization
    parameter Types.Init initType=Types.Init.InitialValues
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
    Medium.EnthalpyFlowRate H_flow_top[nTopPorts]
      "Enthalpy flow rates from the top ports in to the tank";
    Medium.EnthalpyFlowRate port_b_H_flow_bottom[nPorts]
      "Enthalpy flow rates from the bottom ports in to the tank";
    Medium.MassFlowRate mXi_flow_top[nTopPorts, Medium.nXi]
      "Substance mass flow rates from the top ports in to the tank";
    Medium.MassFlowRate port_b_mXi_flowottom[nPorts, Medium.nXi]
      "Substance mass flow rates from the bottom ports in to the tank";
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

  // Only one connection allowed to a port to avoid unwanted ideal mixing
/*
for i in 1:nTopPorts loop
  assert(cardinality(topPorts[i]) <= 1,"
topPorts[" + String(i) + "] of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller.
");
end for;
 
for i in 1:nPorts loop
  assert(cardinality(ports[i]) <= 1,"
ports[" + String(i) + "] of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller.
");
end for;
*/

  // Total quantities
    medium.p = p_ambient;
    V = crossArea*level + V0 "Volume of fluid";
    m = V*medium.d "Mass of fluid";
    mXi = m*medium.Xi "Mass of fluid components";
    U = m*medium.u "Internal energy of fluid";

  // Mass balances
  if dynamicsType < Types.Dynamics.SteadyStateMass then
    der(m) = sum(topPorts.m_flow) + sum(ports.m_flow);
    for i in 1:Medium.nXi loop
      der(mXi[i]) = sum(mXi_flow_top[:,i]) + sum(port_b_mXi_flowottom[:,i]);
    end for;
  else
    0 = sum(topPorts.m_flow) + sum(ports.m_flow);
    for i in 1:Medium.nXi loop
      0 = sum(mXi_flow_top[:,i]) + sum(port_b_mXi_flowottom[:,i]);
    end for;
  end if;

  // Energy balance
  if dynamicsType < Types.Dynamics.SteadyState then
    if Medium.singleState then
      der(U) = sum(H_flow_top) + sum(port_b_H_flow_bottom)
          "Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)";
    else
      der(U) = sum(H_flow_top) + sum(port_b_H_flow_bottom) - p_ambient*der(V);
    end if;
  else
    if Medium.singleState then
      0 = sum(H_flow_top) + sum(port_b_H_flow_bottom)
          "Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)";
    else
      0 = sum(H_flow_top) + sum(port_b_H_flow_bottom) - p_ambient*der(V);
    end if;
  end if;

  // Properties at top ports
    for i in 1:nTopPorts loop
       // It is assumed that fluid flows only from one of the top ports in to the tank and never vice versa
       H_flow_top[i]     = topPorts[i].m_flow*actualStream(topPorts[i].h_outflow);
       mXi_flow_top[i,:] = topPorts[i].m_flow*actualStream(topPorts[i].Xi_outflow);
       topPorts[i].p     = p_ambient;
       topPorts[i].h_outflow = h_start;
/*
       assert(topPorts[i].m_flow > -1, "Mass flows out of tank via topPorts[" + String(i) + "]\n" +
                                         "This indicates a wrong model");
*/
    end for;

  // Properties at bottom ports
    for i in 1:nPorts loop
       port_b_H_flow_bottom[i]   = ports[i].m_flow*actualStream(ports[i].h_outflow);
       port_b_mXi_flowottom[i,:] = ports[i].m_flow*actualStream(ports[i].Xi_outflow);
       aboveLevel[i] = level >= (portsData[i].portLevel + ports_emptyPipeHysteresis[i])
                       or pre(aboveLevel[i]) and level >= (portsData[i].portLevel - ports_emptyPipeHysteresis[i]);
       levelAbovePort[i] = if aboveLevel[i] then level - portsData[i].portLevel else 0;
       ports[i].h_outflow = medium.h;
       ports[i].Xi_outflow = medium.Xi;

       if stiffCharacteristicForEmptyPort then
          // If port is above fluid level, use large zeta if fluid flows out of port (= small mass flow rate)
          zeta_out[i] = 1 + (if aboveLevel[i] then 0 else zetaLarge);
          ports[i].p = p_ambient + levelAbovePort[i]*system.g*medium.d
                               + Modelica_Fluid.Utilities.regSquare2(ports[i].m_flow, m_flow_small,
                                     lossConstant_D_zeta(portsData[i].diameter, 0.01)/medium.d,
                                     lossConstant_D_zeta(portsData[i].diameter, zeta_out[i])/medium.d);
          ports_m_flow_out[i] = false;

       else
          // Handling according to Remelhe/Poschlad
          ports_m_flow_out[i] = (pre(ports_m_flow_out[i]) and not ports[i].p>p_ambient)
                                     or ports[i].m_flow < -1e-6;
         if aboveLevel[i] then
             ports[i].p = p_ambient + levelAbovePort[i]*system.g*medium.d -
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
      Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1},
          initialScale=0.2), graphics={
          Rectangle(
            extent={{-100,-100},{100,100}},
            lineColor={255,255,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent=DynamicSelect({{-100,-100},{100,0}}, {{-100,-100},{100,(-100
                 + 200*level/levelMax)}}),
            lineColor={0,127,255},
            fillColor={85,170,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-94,19},{96,-1}},
            lineColor={0,0,0},
            textString=DynamicSelect(" ", realString(
                level,
                1,
                3))),
          Line(
            points={{-100,100},{100,100}},
            color={0,0,0},
            pattern=LinePattern.Dot),
          Text(
            extent={{-94,90},{95,60}},
            lineColor={0,0,255},
            textString="%name"),
          Text(
            extent={{-95,-85},{95,-65}},
            lineColor={0,0,0},
            textString="%level_start"),
          Text(
            extent={{-95,-55},{95,-35}},
            lineColor={0,0,0},
            textString="level_start ="),
          Text(
            extent={{-95,50},{95,30}},
            lineColor={0,0,0},
            textString="level ="),
          Line(points={{-100,100},{-100,-100},{100,-100},{100,100}}, color={0,0,
                0})}),
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
      Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1},
          initialScale=0.2), graphics),
      uses(Modelica(version="2.2.1"), Modelica_Fluid(version="0.952")));
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
      "Lumped volume with dynamic mass and energy balance"
      import Modelica_Fluid.Types;
        outer Modelica_Fluid.System system "System properties";
        replaceable package Medium = 
          Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
            annotation (choicesAllMatching = true);
        parameter Modelica_Fluid.Types.Dynamics dynamicsType=system.dynamicsType
        "Dynamics option" 
          annotation(Evaluate=true, Dialog(tab = "Assumptions"));
        parameter Types.Init initType=
                  system.initType "Initialization option" 
          annotation(Evaluate=true, Dialog(tab = "Initialization"));
        parameter Medium.AbsolutePressure p_start = system.p_start
        "Start value of pressure" 
          annotation(Dialog(tab = "Initialization"));
        parameter Boolean use_T_start = true
        "= true, use T_start, otherwise h_start" 
          annotation(Dialog(tab = "Initialization"), Evaluate=true);
        parameter Medium.Temperature T_start=
          if use_T_start then system.T_start else Medium.temperature_phX(p_start,h_start,X_start)
        "Start value of temperature" 
          annotation(Dialog(tab = "Initialization", enable = use_T_start));
        parameter Medium.SpecificEnthalpy h_start=
          if use_T_start then Medium.specificEnthalpy_pTX(p_start, T_start, X_start) else Medium.h_default
        "Start value of specific enthalpy" 
          annotation(Dialog(tab = "Initialization", enable = not use_T_start));
        parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default
        "Start value of mass fractions m_i/m" 
          annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));

        Medium.BaseProperties medium(
          preferredMediumStates=true,
          p(start=p_start),
          h(start=h_start),
          T(start=T_start),
          Xi(start=X_start[1:Medium.nXi]));
        SI.Energy U "Internal energy of fluid";
        SI.Mass m "Mass of fluid";
        SI.Mass[Medium.nXi] mXi "Masses of independent components in the fluid";

        // variables that need to be defined by an extending class
        // make Qs_flow an input, allowing its definition and modification with binding equations,
        // e.g. to add a HeatPort to an existing model.
        SI.Volume fluidVolume "Volume";
        SI.MassFlowRate ms_flow "Mass flows across boundaries";
        SI.MassFlowRate[Medium.nXi] msXi_flow
        "Substance mass flows across boundaries";
        SI.EnthalpyFlowRate Hs_flow
        "Enthalpy flow across boundaries or energy source/sink";
        input SI.HeatFlowRate Qs_flow
        "Heat flow across boundaries or energy source/sink";
        SI.Power Ws_flow "Work flow across boundaries or source term";
      equation

        // Total quantities
        m = fluidVolume*medium.d;
        mXi = m*medium.Xi;
        U = m*medium.u;

        // Mass and energy balances
        if dynamicsType < Types.Dynamics.SteadyStateMass then
          der(m) = ms_flow;
          der(mXi) = msXi_flow;
        else
          0 = ms_flow;
          zeros(Medium.nXi) = msXi_flow;
        end if;
        if dynamicsType < Types.Dynamics.SteadyState then
          der(U) = Hs_flow + Qs_flow + Ws_flow;
        else
          0 = Hs_flow + Qs_flow + Ws_flow;
        end if;

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
          if use_T_start then
            der(medium.T) = 0;
          else
            der(medium.h) = 0;
          end if;
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
        annotation (
          Documentation(info="<html>
Base class for an ideally mixed fluid volume with the ability to store mass and energy. 
The following source terms are part of the energy balance and must be specified in an extending class:
<ul>
<li><tt><b>input Qs_flow</b></tt>, e.g. convective or latent heat flow rate across segment boundary, and</li> 
<li><tt><b>Ws_flow</b></tt>, work term, e.g. p*der(fluidVolume) if the volume is not constant.</li>
</ul>
The component volume <tt><b>fluidVolume</b></tt> is a variable which needs to be set in the extending class to complete the model.
<p>
Further source terms must be defined by an extending class for fluid flow across the segment boundary:
</p>
<ul>
<li><tt><b>Hs_flow</b></tt>, enthalpy flow,</li> 
<li><tt><b>ms_flow</b></tt>, mass flow, and</li> 
<li><tt><b>msXi_flow</b></tt>, substance mass flow.</li> 
</ul>
<b>Note:</b>
<p>
<tt>Qs_flow</tt> is defined as input allowing its definition and modification by extending classes, e.g. to add a HeatPort to an existing model.
</p>
<p>
<tt>Hs_flow</tt> is not defined as input as it needs to be defined carefully together with <tt>ms_flow</tt> and modifications can easily lead to bad models. 
Moreover an input might mislead a tool to break equation systems, resulting in inefficient models for the treatment of flow reversal. 
</p>
</html>"),Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
                -100},{100,100}}),
                  graphics));
      end PartialLumpedVolume;

      partial model PartialLumpedVolumePorts
      "Closed volume with a vector of fluid ports"
        extends PartialLumpedVolume;

      //Port definitions
        parameter Integer nPorts(min=1)=1 "Number of ports" annotation(Dialog(group="Ports"));
        parameter SI.Diameter portDiameters[nPorts] = ones(nPorts)
        "Inner (hydraulic) diameters of ports (array)"   annotation(Dialog(group="Ports",enable= not neglectPortDiameters));
        Interfaces.FluidPorts_b[nPorts] ports(
                                      redeclare each package Medium = Medium)
        "Fluid outlets" 
          annotation (Placement(transformation(extent={{-10,-40},{10,40}},
            rotation=-90,
            origin={0,-100}),
            iconTransformation(extent={{-10,40},{10,-40}},
            rotation=-90,
            origin={0,-100})));

      //Transformation of kinetic energy
          parameter Boolean neglectPortDiameters=true
        "=true, kinetic energy and dissipation is accounted for in port pressure"
                                                                                                            annotation(Evaluate=true, Dialog(tab="Assumptions"));
          parameter Real[nPorts] zeta_in=fill(0, nPorts)
        "Hydraulic resistance into volume, 1 for total dissipation of kinetic energy and uniform flow distribution in pipe"
                                                                                                            annotation(Dialog(tab="Assumptions",enable= not neglectPortDiameters));
          parameter Real[nPorts] zeta_out=fill(1, nPorts)
        "Hydraulic resistance out of volume, 0 for ideal smooth outlet"   annotation(Dialog(tab="Assumptions",enable= not neglectPortDiameters));

        Medium.EnthalpyFlowRate ports_H_flow[nPorts];
        Medium.MassFlowRate ports_mXi_flow[nPorts,Medium.nXi];
        Medium.MassFlowRate[Medium.nXi] sum_ports_mXi_flow
        "Substance mass flows through ports";
        Medium.ExtraPropertyFlowRate ports_mC_flow[nPorts,Medium.nC];

    protected
        parameter SI.Area[nPorts] portArea=Modelica.Constants.pi/4*{portDiameters[i]^2 for i in 1:nPorts};

      equation
        ms_flow = sum(ports.m_flow);
        msXi_flow = sum_ports_mXi_flow;
        Hs_flow = sum(ports_H_flow);

        // Only one connection allowed to a port to avoid unwanted ideal mixing
        for i in 1:nPorts loop
          assert(cardinality(ports[i]) <= 1,"
each ports[i] of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections, which is usually not the intention
of the modeller. Increase nPorts to add an additional port.
");
        end for;

        // Boundary conditions
        for i in 1:nPorts loop
          if neglectPortDiameters then
            ports[i].p = medium.p;
          else
            ports[i].p = medium.p - smooth(2, noEvent(ports[i].m_flow^2/(2*medium.d*
                  portArea[i]^2)*(if ports[i].m_flow < 0 then (1 + zeta_out[i]) else (1
                   - zeta_in[i]))));
          end if;
        end for;
        ports.h_outflow = fill(medium.h, nPorts);
        ports.Xi_outflow = fill(medium.Xi, nPorts);

        for i in 1:nPorts loop
          ports_H_flow[i] = ports[i].m_flow * actualStream(ports[i].h_outflow)
          "Enthalpy flow";
          ports_mXi_flow[i,:] = ports[i].m_flow * actualStream(ports[i].Xi_outflow)
          "Component mass flow";
          ports_mC_flow[i,:] = ports[i].m_flow * actualStream(ports[i].C_outflow)
          "Trace substance mass flow";
        end for;
        for i in 1:Medium.nXi loop
          sum_ports_mXi_flow[i] = sum(ports_mXi_flow[:,i]);
        end for;
        // Steady-state trace substance mass balance
        for i in 1:Medium.nC loop
          sum(ports_mC_flow[:,i]) = 0;
        end for;

       annotation (
        Documentation(info="<html>
This base class extends PartialLumpedVolume by adding a vector of fluid ports 
and defining the respective source terms
<ul>
<li><tt>Hs_flow</tt>, enthalpy flow,</li> 
<li><tt>ms_flow</tt>, mass flow, and</li> 
<li><tt>msXi_flow</tt>, substance mass flow.</li> 
</ul>
An extending class still needs to define:
<ul>
<li><tt>Qs_flow</tt>, e.g. convective or latent heat flow rate across segment boundary,</li> 
<li><tt>Ws_flow</tt>, work term, e.g. p*der(V) if the volume is not constant, and</li>
<li><tt>V_lumped</tt>, the volume of the segment.</li>
</ul>
</html>"),
        Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},
                {100,100}}),
                graphics),
          Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},
                {100,100}}), graphics={Text(
              extent={{-150,110},{150,150}},
              textString="%name",
              lineColor={0,0,255}), Text(
              extent={{-150,-110},{150,-140}},
              lineColor={0,0,0},
              textString="V=%V")}));
      end PartialLumpedVolumePorts;

  partial model PartialDistributedVolume "Base class for a finite volume model"
    import Modelica_Fluid.Types;
    outer Modelica_Fluid.System system "System properties";

    replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching = true);

    // Discretization
    parameter Integer n=1 "Number of discrete flow volumes";

    // Assumptions
    parameter Modelica_Fluid.Types.Dynamics dynamicsType=system.dynamicsType
        "Dynamics option" 
      annotation(Evaluate=true, Dialog(tab = "Assumptions"));

    //final parameter Boolean static = dynamicsType == Types.Dynamics.SteadyState
    //  "= true, static balances, no mass or energy is stored" annotation 2;

    //Initialization
    parameter Types.Init initType=system.initType "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));

    parameter Medium.AbsolutePressure p_a_start=system.p_start
        "Start value of pressure at port a" 
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_b_start=p_a_start
        "Start value of pressure at port b" 
      annotation(Dialog(tab = "Initialization"));
    final parameter Medium.AbsolutePressure[n] p_start=if n > 1 then linspace(
          p_a_start - (p_a_start - p_b_start)/(2*n),
          p_b_start + (p_a_start - p_b_start)/(2*n),
          n) else {(p_a_start + p_b_start)/2} "Start value of pressure";

    parameter Boolean use_T_start=true "Use T_start if true, otherwise h_start"
       annotation(Evaluate=true, Dialog(tab = "Initialization"));
    parameter Medium.Temperature T_start=if use_T_start then system.T_start else 
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

    // Total quantities
    SI.Energy[n] U "Internal energy of fluid";
    SI.Mass[n] m "Fluid mass";
    SI.Mass[n,Medium.nXi] mXi "Substance mass";

    Medium.BaseProperties[n] medium(
      each preferredMediumStates=if (dynamicsType == Types.Dynamics.SteadyState) then false else true,
      p(start=p_start),
      each h(start=h_start),
      each T(start=T_start),
      each Xi(start=X_start[1:Medium.nXi]));

    //Source terms, have to be set in inheriting class (to zero if not used)
    SI.Volume[n] fluidVolume
        "Discretized volume, determine in inheriting class";
    Medium.MassFlowRate[n] ms_flow "Mass flow rate, source or sink";
    Medium.MassFlowRate[n,Medium.nXi] msXi_flow
        "Independent mass flow rates, source or sink";
    SI.EnthalpyFlowRate[n] Hs_flow "Enthalpy flow rate, source or sink";
    input SI.HeatFlowRate[n] Qs_flow "Heat flow rate, source or sink";
    SI.Power[n] Ws_flow "Mechanical power, p*der(V) etc.";

  equation
    // Total quantities
    for i in 1:n loop
      m[i] =fluidVolume[i]*medium[i].d;
      mXi[i, :] = m[i]*medium[i].Xi;
      U[i] = m[i]*medium[i].u;
    end for;

    // Mass and energy balances
    if dynamicsType < Types.Dynamics.SteadyStateMass then
    // dynamic mass balances, n "thermal" states, n pressure states (if not singleState_hydraulic)
      for i in 1:n loop
        der(m[i]) = ms_flow[i];
        der(mXi[i, :]) = msXi_flow[i, :];
      end for;
    else
    // steady state mass balances, no numerical states, no flow reversal possible
      for i in 1:n loop
        0 = ms_flow[i];
        zeros(Medium.nXi) = msXi_flow[i, :];
      end for;
    end if;
    if dynamicsType < Types.Dynamics.SteadyState then
    // dynamic energy balances, n "thermal" states, n pressure states (if not singleState_hydraulic)
      for i in 1:n loop
        der(U[i]) = Hs_flow[i] + Ws_flow[i] + Qs_flow[i];
      end for;
    else
    // steady state energy balances, no numerical states, no flow reversal possible
      for i in 1:n loop
        0 = Hs_flow[i] + Ws_flow[i] + Qs_flow[i];
      end for;
    end if;

  initial equation
    // Initial conditions
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

     annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
                -100},{100,100}}),
                         graphics),
                          Icon(coordinateSystem(preserveAspectRatio=true,
              extent={{-100,-100},{100,100}}), graphics),
        Documentation(info="<html>
Base class for <tt><b>n</b></tt> ideally mixed fluid volumes with the ability to store mass and energy.
It is inteded to model a one-dimensional spatial discretization of fluid flow according to the finite volume method. 
The following source terms are part of the energy balance and must be specified in an extending class:
<ul>
<li><tt><b>input Qs_flow[n]</b></tt>, e.g. convective or latent heat flow rate across segment boundary, and</li> 
<li><tt><b>Ws_flow[n]</b></tt>, work term, e.g. p*der(fluidVolume) if the volume is not constant.</li>
</ul>
The component volume <tt><b>fluidVolume[n]</b></tt> is a variable which needs to be set in the extending class to complete the model.
<p>
Further source terms must be defined by an extending class for fluid flow across the segment boundary:
</p>
<ul>
<li><tt><b>Hs_flow[n]</b></tt>, enthalpy flow,</li> 
<li><tt><b>ms_flow[n]</b></tt>, mass flow, and</li> 
<li><tt><b>msXi_flow[n]</b></tt>, substance mass flow.</li> 
</ul>
<b>Note:</b>
<p>
<tt>Qs_flow</tt> is defined as input allowing its definition and modification by extending classes, e.g. to add a HeatPort to an existing model.
</p>
<p>
<tt>Hs_flow</tt> is not defined as input as it needs to be defined carefully together with <tt>ms_flow</tt> and modifications can easily lead to bad models. 
Moreover an input might mislead a tool to break equation systems, resulting in inefficient models for the treatment of flow reversal. 
</p>
</html>"));
  end PartialDistributedVolume;
  end BaseClasses;
  annotation (Documentation(info="<html>
 
</html>"));
end Volumes;
