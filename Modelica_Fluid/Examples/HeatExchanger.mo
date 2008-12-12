within Modelica_Fluid.Examples;
package HeatExchanger "Demo of a heat exchanger model"
  model HeatExchangerSimulation "simulation for the heat exchanger model"

  extends Modelica.Icons.Example;

  //replaceable package Medium = Modelica.Media.Water.StandardWater;
  package Medium = Modelica.Media.Incompressible.Examples.Essotherm650;
    Modelica_Fluid.Examples.HeatExchanger.BaseClasses.BasicHX HEX(
      c_wall=500,
      use_T_start=true,
      nNodes=20,
      length=2,
      m_flow_start_1=0.2,
      m_flow_start_2=0.2,
      redeclare model PressureLoss_1 = 
          Modelica_Fluid.Pipes.BaseClasses.PressureLoss.DetailedFlow (
              use_d_nominal=true,use_eta_nominal=true,eta_nominal=0.01),
      redeclare model PressureLoss_2 = 
          Modelica_Fluid.Pipes.BaseClasses.PressureLoss.DetailedFlow (
              use_d_nominal=true,use_eta_nominal=true,eta_nominal=0.01),
      k_wall=100,
      initType=Modelica_Fluid.Types.Init.SteadyStateHydraulic,
      s_wall=0.005,
      crossArea_1=4.5e-4,
      crossArea_2=4.5e-4,
      perimeter_1=0.075,
      perimeter_2=0.075,
      area_h_1=0.075*2*20,
      area_h_2=0.075*2*20,
      d_wall=900,
      redeclare package Medium_1 = 
          Medium,
      redeclare package Medium_2 = 
          Medium,
      redeclare model HeatTransfer_1 = 
          Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_constAlpha (alpha0=
             1000),
      Twall_start=300,
      dT=10,
      T_start_1=304,
      T_start_2=300)       annotation (Placement(transformation(extent={{
              -26,-14},{34,46}}, rotation=0)));

    Modelica_Fluid.Sources.FixedBoundary_pTX ambient2(
      p=1e5,
      T=280,
      redeclare package Medium = Medium)                              annotation (Placement(
          transformation(extent={{82,-28},{62,-8}}, rotation=0)));
    Modelica_Fluid.Sources.FixedBoundary_pTX ambient1(
      p=1e5,
      T=300,
      redeclare package Medium = Medium)                              annotation (Placement(
          transformation(extent={{82,24},{62,44}}, rotation=0)));
    Modelica_Fluid.Sources.PrescribedMassFlowRate_TX massFlowRate2(
      m_flow=0.2,
      T=360,
      redeclare package Medium = Medium,
      useFlowRateInput=true,
      useTemperatureInput=false,
      useCompositionInput=false) 
                  annotation (Placement(transformation(extent={{-66,24},{-46,44}},
            rotation=0)));
    Modelica_Fluid.Sources.PrescribedMassFlowRate_TX massFlowRate1(
      T=300,
      m_flow=0.5,
      redeclare package Medium = Medium) 
                   annotation (Placement(transformation(extent={{-66,-10},{-46,10}},
            rotation=0)));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
              -100},{100,100}}),
                        graphics),
                         experiment(StopTime=100, Tolerance=1e-005));
    Modelica.Blocks.Sources.Ramp Ramp1(
      startTime=50,
      duration=5,
      height=-1,
      offset=0.5)   annotation (Placement(transformation(extent={{-100,24},{-80,
              44}}, rotation=0)));
    inner Modelica_Fluid.System system 
                                     annotation (Placement(transformation(extent=
              {{60,70},{80,90}}, rotation=0)));
  equation
    connect(massFlowRate1.ports[1], HEX.port_a1)        annotation (Line(points={
            {-46,0},{-40,0},{-40,15.4},{-29,15.4}}, color={0,127,255}));
    connect(HEX.port_b1, ambient1.ports[1])        annotation (Line(points={{37,
            15.4},{48.5,15.4},{48.5,34},{62,34}}, color={0,127,255}));
    connect(Ramp1.y, massFlowRate2.m_flow_in) annotation (Line(points={{-79,34},{
            -74,34},{-74,42},{-66,42}},   color={0,0,127}));
    connect(massFlowRate2.ports[1], HEX.port_b2) 
                                             annotation (Line(
        points={{-46,34},{-40,34},{-40,29.8},{-29,29.8}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(HEX.port_a2, ambient2.ports[1]) 
                                        annotation (Line(
        points={{37,2.2},{42,2},{50,2},{50,-18},{62,-18}},
        color={0,127,255},
        smooth=Smooth.None));
  end HeatExchangerSimulation;

  package BaseClasses "Additional models for heat exchangers"
    extends Modelica_Fluid.Icons.BaseClassLibrary;

    model BasicHX "Simple heat exchanger model"
      outer Modelica_Fluid.System system "System properties";
      //General
      parameter Integer nNodes(min=1) = 1 "Spatial segmentation";
      replaceable package Medium_1 = Modelica.Media.Water.StandardWater constrainedby
        Modelica.Media.Interfaces.PartialMedium "Fluid 1" 
                                                        annotation(choicesAllMatching, Dialog(tab="General",group="Fluid 1"));
      replaceable package Medium_2 = Modelica.Media.Water.StandardWater constrainedby
        Modelica.Media.Interfaces.PartialMedium "Fluid 2" 
                                                        annotation(choicesAllMatching,Dialog(tab="General", group="Fluid 2"));
      parameter SI.Area crossArea_1 "Cross sectional area" annotation(Dialog(tab="General",group="Fluid 1"));
      parameter SI.Area crossArea_2 "Cross sectional area" annotation(Dialog(tab="General",group="Fluid 2"));
      parameter SI.Length perimeter_1 "Flow channel perimeter" annotation(Dialog(tab="General",group="Fluid 1"));
      parameter SI.Length perimeter_2 "Flow channel perimeter" annotation(Dialog(tab="General",group="Fluid 2"));
      parameter SI.Length length(min=0) "Length of flow path for both fluids";
      parameter SI.Length s_wall(min=0) "Wall thickness";
      // Heat transfer
      replaceable model HeatTransfer_1 = 
          Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_constAlpha 
        constrainedby
        Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialFlowHeatTransfer
        "Heat transfer model" annotation(choicesAllMatching, Dialog(tab="General", group="Fluid 1"));

      replaceable model HeatTransfer_2 = 
          Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_constAlpha 
        constrainedby
        Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialFlowHeatTransfer
        "Heat transfer model" annotation(choicesAllMatching, Dialog(tab="General", group="Fluid 2"));

      parameter SI.Area area_h_1 "Heat transfer area" annotation(Dialog(tab="General",group="Fluid 1"));
      parameter SI.Area area_h_2 "Heat transfer area" annotation(Dialog(tab="General",group="Fluid 2"));
     //Wall
      parameter SI.Density d_wall "Density of wall material" annotation(Dialog(tab="General", group="Solid material properties"));
      parameter SI.SpecificHeatCapacity c_wall
        "Specific heat capacity of wall material" annotation(Dialog(tab="General", group="Solid material properties"));
      final parameter SI.Area area_h=(area_h_1 + area_h_2)/2
        "Heat transfer area";
      final parameter SI.Mass m_wall=d_wall*area_h*s_wall "Wall mass";
      parameter SI.ThermalConductivity k_wall
        "Thermal conductivity of wall material" 
        annotation (Dialog(group="Solid material properties"));

      // Assumptions
      parameter Boolean allowFlowReversal = system.allowFlowReversal
        "allow flow reversal, false restricts to design direction (port_a -> port_b)"
        annotation(Dialog(tab="Assumptions"), Evaluate=true);
      parameter Modelica_Fluid.Types.Dynamics dynamicsType=system.dynamicsType
        "Dynamics option" 
        annotation(Evaluate=true, Dialog(tab = "Assumptions"));

      //Initialization pipe 1
      parameter Types.Init initType=Types.Init.InitialValues
        "Initialization option" 
        annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter SI.Temperature Twall_start "Start value of wall temperature" 
                                                                            annotation(Dialog(tab="Initialization", group="Wall"));
      parameter SI.Temperature dT "Start value for pipe_1.T - pipe_2.T" 
        annotation (Dialog(tab="Initialization", group="Wall"));
      parameter Boolean use_T_start=true
        "Use T_start if true, otherwise h_start" 
        annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Medium_1.AbsolutePressure p_a_start1=Medium_1.p_default
        "Start value of pressure" 
        annotation(Dialog(tab = "Initialization", group = "Fluid 1"));
      parameter Medium_1.AbsolutePressure p_b_start1=Medium_1.p_default
        "Start value of pressure" 
        annotation(Dialog(tab = "Initialization", group = "Fluid 1"));
      parameter Medium_1.Temperature T_start_1=if use_T_start then Medium_1.
          T_default else Medium_1.temperature_phX(
            (p_a_start1 + p_b_start1)/2,
            h_start_1,
            X_start_1) "Start value of temperature" 
        annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Fluid 1", enable = use_T_start));
      parameter Medium_1.SpecificEnthalpy h_start_1=if use_T_start then Medium_1.specificEnthalpy_pTX(
            (p_a_start1 + p_b_start1)/2,
            T_start_1,
            X_start_1[1:Medium_1.nXi]) else Medium_1.h_default
        "Start value of specific enthalpy" 
        annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Fluid 1", enable = not use_T_start));
      parameter Medium_1.MassFraction X_start_1[Medium_1.nX]=Medium_1.X_default
        "Start value of mass fractions m_i/m" 
        annotation (Dialog(tab="Initialization", group = "Fluid 1", enable=(Medium_1.nXi > 0)));
      parameter Medium_1.MassFlowRate m_flow_start_1 = system.m_flow_start
        "Start value of mass flow rate" annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Fluid 1"));
      //Initialization pipe 2

      parameter Medium_2.AbsolutePressure p_a_start2=Medium_2.p_default
        "Start value of pressure" 
        annotation(Dialog(tab = "Initialization", group = "Fluid 2"));
      parameter Medium_2.AbsolutePressure p_b_start2=Medium_2.p_default
        "Start value of pressure" 
        annotation(Dialog(tab = "Initialization", group = "Fluid 2"));
      parameter Medium_2.Temperature T_start_2=if use_T_start then Medium_2.
          T_default else Medium_2.temperature_phX(
            (p_a_start2 + p_b_start2)/2,
            h_start_2,
            X_start_2) "Start value of temperature" 
        annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Fluid 2", enable = use_T_start));
      parameter Medium_2.SpecificEnthalpy h_start_2=if use_T_start then Medium_2.specificEnthalpy_pTX(
            (p_a_start2 + p_b_start2)/2,
            T_start_2,
            X_start_2[1:Medium_2.nXi]) else Medium_2.h_default
        "Start value of specific enthalpy" 
        annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Fluid 2", enable = not use_T_start));
      parameter Medium_2.MassFraction X_start_2[Medium_2.nX]=Medium_2.X_default
        "Start value of mass fractions m_i/m" 
        annotation (Dialog(tab="Initialization", group = "Fluid 2", enable=Medium_2.nXi>0));
      parameter Medium_2.MassFlowRate m_flow_start_2 = system.m_flow_start
        "Start value of mass flow rate"    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Fluid 2"));

      //Pressure drop and heat transfer
      replaceable model PressureLoss_1 = 
          Modelica_Fluid.Pipes.BaseClasses.PressureLoss.QuadraticTurbulentFlow 
        constrainedby
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.PartialFlowPressureLoss
        "Characteristic of wall friction"                                                                                                   annotation(choicesAllMatching, Dialog(tab="General", group="Fluid 1"));
      replaceable model PressureLoss_2 = 
          Modelica_Fluid.Pipes.BaseClasses.PressureLoss.QuadraticTurbulentFlow 
        constrainedby
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.PartialFlowPressureLoss
        "Characteristic of wall friction"                                                                                                   annotation(choicesAllMatching, Dialog(tab="General", group="Fluid 2"));
      parameter SI.Length roughness_1=2.5e-5
        "Absolute roughness of pipe (default = smooth steel pipe)" annotation(Dialog(tab="General", group="Fluid 1"));
      parameter SI.Length roughness_2=2.5e-5
        "Absolute roughness of pipe (default = smooth steel pipe)" annotation(Dialog(tab="General", group="Fluid 2"));

      //Display variables
      SI.HeatFlowRate Q_flow_1 "Total heat flow rate of pipe 1";
      SI.HeatFlowRate Q_flow_2 "Total heat flow rate of pipe 2";

      BaseClasses.WallConstProps wall(
        d_wall=d_wall,
        c_wall=c_wall,
        T_start=Twall_start,
        k_wall=k_wall,
        dT=dT,
        s=s_wall,
        initType=initType,
        n=nNodes,
        area_h=(crossArea_1 + crossArea_2)/2) 
        annotation (Placement(transformation(extent={{-29,-23},{9,35}},  rotation=
               0)));

      Modelica_Fluid.Pipes.DistributedPipe pipe_1(
        redeclare package Medium = Medium_1,
        isCircular=false,
        diameter=0,
        nNodes=nNodes,
        allowFlowReversal=allowFlowReversal,
        dynamicsType=dynamicsType,
        length=length,
        redeclare model HeatTransfer = HeatTransfer_1(area=area_h_1),
        initType=initType,
        use_T_start=use_T_start,
        T_start=T_start_1,
        h_start=h_start_1,
        X_start=X_start_1,
        m_flow_start=m_flow_start_1,
        perimeter=perimeter_1,
        crossArea=crossArea_1,
        roughness=roughness_1,
        redeclare model PressureLoss = PressureLoss_1)   annotation (Placement(transformation(extent={{-40,-80},
                {20,-20}},        rotation=0)));

      Modelica_Fluid.Pipes.DistributedPipe pipe_2(
        redeclare package Medium = Medium_2,
        nNodes=nNodes,
        allowFlowReversal=allowFlowReversal,
        dynamicsType=dynamicsType,
        length=length,
        isCircular=false,
        diameter=0,
        redeclare model HeatTransfer = HeatTransfer_2(area=area_h_2),
        use_T_start=use_T_start,
        T_start=T_start_2,
        h_start=h_start_2,
        X_start=X_start_2,
        initType=initType,
        m_flow_start=m_flow_start_2,
        perimeter=perimeter_2,
        crossArea=crossArea_2,
        p_a_start=p_a_start1,
        p_b_start=p_b_start2,
        roughness=roughness_2,
        redeclare model PressureLoss = PressureLoss_2) 
                  annotation (Placement(transformation(extent={{20,88},{-40,28}},
              rotation=0)));

      annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
                -100},{100,100}},
            grid={1,1}),  graphics),
                           Icon(coordinateSystem(preserveAspectRatio=false,
              extent={{-100,-100},{100,100}},
            grid={1,1}), graphics={
            Rectangle(
              extent={{-100,-26},{100,-30}},
              lineColor={0,0,0},
              fillColor={95,95,95},
              fillPattern=FillPattern.Forward),
            Rectangle(
              extent={{-100,30},{100,26}},
              lineColor={0,0,0},
              fillColor={95,95,95},
              fillPattern=FillPattern.Forward),
            Rectangle(
              extent={{-100,60},{100,30}},
              lineColor={0,0,0},
              fillPattern=FillPattern.HorizontalCylinder,
              fillColor={0,63,125}),
            Rectangle(
              extent={{-100,-30},{100,-60}},
              lineColor={0,0,0},
              fillPattern=FillPattern.HorizontalCylinder,
              fillColor={0,63,125}),
            Rectangle(
              extent={{-100,26},{100,-26}},
              lineColor={0,0,0},
              fillPattern=FillPattern.HorizontalCylinder,
              fillColor={0,128,255}),
            Text(
              extent={{-150,110},{150,70}},
              lineColor={0,0,255},
              textString="%name"),
            Line(
              points={{30,-85},{-60,-85}},
              color={0,128,255},
              smooth=Smooth.None),
            Polygon(
              points={{20,-70},{60,-85},{20,-100},{20,-70}},
              lineColor={0,128,255},
              smooth=Smooth.None,
              fillColor={0,128,255},
              fillPattern=FillPattern.Solid),
            Line(
              points={{30,77},{-60,77}},
              color={0,128,255},
              smooth=Smooth.None),
            Polygon(
              points={{-50,92},{-90,77},{-50,62},{-50,92}},
              lineColor={0,128,255},
              smooth=Smooth.None,
              fillColor={0,128,255},
              fillPattern=FillPattern.Solid)}),
        Documentation(info="<html>
Simple model of a heat exchanger consisting of two pipes and one wall in between. 
For both fluids geometry parameters, such as heat transfer area and cross section as well as heat transfer and pressure drop correlations may be chosen. 
The flow scheme may be concurrent or counterflow, defined by the respective flow directions of the fluids entering the component.
The design flow direction with positive m_flow variables is counterflow.
</html>"));
      Modelica_Fluid.Interfaces.FluidPort_b port_b1(redeclare package Medium = 
            Medium_1) annotation (Placement(transformation(extent={{100,-12},{120,
                8}}, rotation=0)));
      Modelica_Fluid.Interfaces.FluidPort_a port_a1(redeclare package Medium = 
            Medium_1) annotation (Placement(transformation(extent={{-120,-12},{
                -100,8}}, rotation=0)));
      Modelica_Fluid.Interfaces.FluidPort_b port_b2(redeclare package Medium = 
            Medium_2) annotation (Placement(transformation(extent={{-120,36},{
                -100,56}}, rotation=0)));
      Modelica_Fluid.Interfaces.FluidPort_a port_a2(redeclare package Medium = 
            Medium_2) annotation (Placement(transformation(extent={{100,-56},{120,
                -36}}, rotation=0)));

    equation
      Q_flow_1 = sum(pipe_1.heatTransfer.Q_flow);
      Q_flow_2 = sum(pipe_2.heatTransfer.Q_flow);
      connect(pipe_2.port_b, port_b2) annotation (Line(
          points={{-40,58},{-76,58},{-76,46},{-110,46}},
          color={0,127,255},
          thickness=0.5));
      connect(pipe_1.port_b, port_b1) annotation (Line(
          points={{20,-50},{42,-50},{42,-2},{110,-2}},
          color={0,127,255},
          thickness=0.5));
      connect(pipe_1.port_a, port_a1) annotation (Line(
          points={{-40,-50},{-75.3,-50},{-75.3,-2},{-110,-2}},
          color={0,127,255},
          thickness=0.5));
      connect(pipe_2.port_a, port_a2) annotation (Line(
          points={{20,58},{65,58},{65,-46},{110,-46}},
          color={0,127,255},
          thickness=0.5));
      connect(wall.heatPort_b, pipe_1.heatPorts) annotation (Line(
          points={{-10,-8.5},{-10,-34.4},{-9.7,-34.4}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(pipe_2.heatPorts[nNodes:-1:1], wall.heatPort_a[1:nNodes]) 
        annotation (Line(
          points={{-10.3,42.4},{-10.3,31.7},{-10,31.7},{-10,20.5}},
          color={127,0,0},
          smooth=Smooth.None));
    end BasicHX;

    model WallConstProps
      "Pipe wall with capacitance, assuming 1D heat conduction and constant material properties"
      parameter Integer n(min=1)=1
        "Segmentation perpendicular to heat conduction";
    //Geometry
      parameter SI.Length s "Wall thickness";
      parameter SI.Area area_h "Heat transfer area";
    //Material properties
      parameter SI.Density d_wall "Density of wall material";
      parameter SI.SpecificHeatCapacity c_wall
        "Specific heat capacity of wall material";
      parameter SI.ThermalConductivity k_wall
        "Thermal conductivity of wall material";
      parameter SI.Mass[n] m=fill(d_wall*area_h*s/n,n)
        "Distribution of wall mass";
    //Initialization
      parameter Types.Init initType=Types.
            Init.NoInit "Initialization option" 
        annotation(Evaluate=true);
      parameter SI.Temperature T_start "Wall temperature start value";
      parameter SI.Temperature dT "Start value for port_b.T - port_a.T";
    //Temperatures
      SI.Temperature[n] Tb(each start=T_start+0.5*dT);
      SI.Temperature[n] Ta(each start=T_start-0.5*dT);
      SI.Temperature[n] T(start=ones(n)*T_start, stateSelect=StateSelect.prefer)
        "Wall temperature";
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] heatPort_a
        "Thermal port" 
        annotation (Placement(transformation(extent={{-20,40},{20,60}}, rotation=0)));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] heatPort_b
        "Thermal port" 
        annotation (Placement(transformation(extent={{-20,-40},{20,-60}}, rotation=
                  0)));

    initial equation
      if initType == Types.Init.SteadyState or initType == Types.Init.SteadyStateHydraulic then
        der(T) = zeros(n);
      elseif initType == Types.Init.InitialValues then
        T = ones(n)*T_start;
      end if;
    equation
      for i in 1:n loop
       assert(m[i]>0, "Wall has negative dimensions");
       c_wall*m[i]*der(T[i]) = heatPort_a[i].Q_flow + heatPort_b[i].Q_flow;
       heatPort_a[i].Q_flow=k_wall/s*(Ta[i]-T[i])*area_h/n;
       heatPort_b[i].Q_flow=k_wall/s*(Tb[i]-T[i])*area_h/n;
      end for;
      Ta=heatPort_a.T;
      Tb=heatPort_b.T;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Rectangle(
              extent={{-100,40},{100,-40}},
              lineColor={0,0,0},
              fillColor={95,95,95},
              fillPattern=FillPattern.Forward), Text(
              extent={{-82,18},{76,-18}},
              lineColor={0,0,0},
              fillColor={0,0,0},
              fillPattern=FillPattern.Forward,
              textString="%name")}),
                                Documentation(revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>",     info="<html>
Simple model of circular (or any other closed shape) wall to be used for pipe (or duct) models. Heat conduction is regarded one dimensional, capacitance is lumped at the arithmetic mean temperature. The spatial discretization (parameter <tt>n</tt>) is meant to correspond to a connected fluid model discretization.
</html>"));
    end WallConstProps;
  end BaseClasses;
end HeatExchanger;