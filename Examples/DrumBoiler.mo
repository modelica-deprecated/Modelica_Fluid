package DrumBoiler 
  "Drum boiler example, see Franke, Rode, Krueger: On-line Optimization of Drum Boiler Startup, 3rd International Modelica Conference, Linkoping, 2003" 
  
  model DrumBoilerSimulation "Simulate start-up of DrumBoiler" 
    extends Modelica.Icons.Example;
    Components.DrumBoiler drumBoiler annotation (extent=[-20, -40; 40, 20]);
    Modelica.Blocks.Sources.TimeTable q_F_Tab(table=[0, 0; 3600, 400; 7210,
          400]) annotation (extent=[-80,0; -60,20]);
    Modelica.Blocks.Sources.TimeTable Y_Valve_Tab(table=[0, 1; 3600, 1; 7210,
           1]) annotation (extent=[-80, -40; -60, -20]);
    annotation (
      Diagram,
      experiment(StopTime=7200),
      Documentation(info="<HTML>
<p>
Apply a ramp to fuel input and hold outlet valve open.
Simulate for 7200 seconds.
</p>
</HTML>"));
  equation 
    connect(q_F_Tab.y,       drumBoiler.q_F) annotation (points=[-59,10; -40,10; 
          -40,-28; -21.35,-28],       style(rgbcolor={0,0,127}));
    connect(Y_Valve_Tab.y,       drumBoiler.Y_Valve) annotation (points=[-59,
          -30; -44,-30; -44,-34; -21.35,-34],    style(
        rgbcolor={0,0,127},
        fillColor=7,
        fillPattern=1));
  end DrumBoilerSimulation;
  
  package Components 
    
    model Evaporator 
      "Simple Evaporator with two states, see Astroem, Bell: Drum-boiler dynamics, Automatica 36, 2000, pp.363-378" 
      
      import Modelica_Fluid.Interfaces.*;
      import Modelica.SIunits.Conversions.*;
      import SI = Modelica.SIunits;
      
      // property and interface declarations
      package Medium = WaterPhaseBoundaryIF97;
      Medium.BaseProperties medium_a(region=1, p=port_a.p) "Medium in port_a";
      Medium.BaseProperties medium_b(region=2, p=port_b.p) "Medium in port_b";
      FluidPort_a port_a(redeclare package Medium = Medium) 
        annotation (extent=[-120, -10; -100, 10]);
      FluidPort_b port_b(redeclare package Medium = Medium) 
        annotation (extent=[120, -10; 100, 10]);
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort 
        annotation (extent=[-10, -120; 10, -100]);
      Modelica.Blocks.Interfaces.RealOutput V(
                                           redeclare type SignalType = 
            SI.Volume) "liquid volume (level)" 
        annotation (extent=[30, 100; 50, 120], rotation=90);
      Modelica.Blocks.Interfaces.RealOutput sigma_D "Thermal stress in metal" 
        annotation (extent=[100, 40; 120, 60]);
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
            extent=[-120, 117; 116, 51],
            string="%name",
            style(gradient=2, fillColor=69)),
          Line(points=[0, -60; 0, -100], style(color=42)),
          Line(points=[40, 99; 40, 60])));
      
      // public parameters
      parameter SI.Mass m_D=300e3 "mass of surrounding drum metal";
      parameter SI.SpecificHeatCapacity cp_D=500 
        "specific heat capacity of drum metal";
      parameter SI.Volume V_t=100 "total volume inside drum";
      parameter SI.Pressure p_start=from_bar(1) "initial pressure";
      parameter SI.Volume V_start=67 "initial liquid volume";
      
    protected 
      SI.Pressure p(start=p_start, stateSelect=StateSelect.prefer) 
        "pressure inside drum boiler";
      SI.Volume V_v "volume of vapour phase";
      SI.Volume V_l(start=V_start, stateSelect=StateSelect.prefer) 
        "volumes of liquid phase";
      SI.SpecificEnthalpy h_v=medium_b.h "specific enthalpy of vapour";
      SI.SpecificEnthalpy h_l=medium_a.h "specific enthalpy of liquid";
      SI.Density rho_v=medium_b.d "density in vapour phase";
      SI.Density rho_l=medium_a.d "density in liquid phase";
      SI.Mass m "total mass of drum boiler";
      SI.Energy U "internal energy";
      SI.Temperature T_D=heatPort.T "temperature of drum";
      SI.HeatFlowRate q_F=heatPort.Q_flow "heat flow rate from furnace";
      SI.SpecificEnthalpy h_W=port_a.h "feed water enthalpy";
      SI.SpecificEnthalpy h_S=medium_b.h "steam enthalpy";
      SI.MassFlowRate qm_W=port_a.m_flow "feed water mass flow rate";
      SI.MassFlowRate qm_S=port_b.m_flow "steam mass flow rate";
    equation 
      
      // balance equations  
      m = rho_v*V_v + rho_l*V_l + m_D;
      U = rho_v*V_v*h_v + rho_l*V_l*h_l - p*V_t + m_D*cp_D*T_D;
      der(m) = qm_W + qm_S;
      der(U) = q_F + qm_W*h_W + qm_S*h_S;
      T_D = medium_a.T;
      // ideal heat transfer between metal and water
      V_t = V_l + V_v;
      
      // pressure and specific total enthalpies at ports
      port_a.p = p;
      port_b.p = p;
      port_b.H_flow = semiLinear(port_b.m_flow, port_b.h, h_v);
      port_a.H_flow = semiLinear(port_a.m_flow, port_a.h, h_l);
      
      // thermal stress
      sigma_D           = -60*der(T_D);
      
      // liquid level 
      V           = V_l;
    end Evaporator;
    
    model Valve "Simple controlled valve with linear pressure drop coefficient" 
      
      import SI = Modelica.SIunits;
      import Modelica.SIunits.Conversions.*;
      import Modelica_Fluid.*;
      extends Interfaces.PartialTwoPortTransport;
      SI.Pressure dp "Pressure loss due to friction";
      Real residue=port_a.p - port_b.p - dp 
        "momentum balance (may be modified)";
      
      parameter Real k=1e-5 "linear valve coefficient";
      
      Modelica.Blocks.Interfaces.RealInput Y "Valve position" 
        annotation (extent=[-10, -80; 10, -60], rotation=90);
      
      annotation (Icon(
          Text(
            extent=[-126, -76; 130, -110],
            style(color=0),
            string=""),
          Text(
            extent=[-120, 130; 116, 64],
            string="%name",
            style(gradient=2, fillColor=69)),
          Line(points=[-60, -50; -60, 50; 60, -50; 60, 50; -60, -50], style(
                color=0, thickness=2)),
          Line(points=[-60, 0; -100, 0], style(color=69)),
          Line(points=[60, 0; 100, 0], style(color=69)),
          Line(points=[0, 0; 0, -72])));
    equation 
      residue = 0;
      port_a.m_flow =Y           *k*dp;
    end Valve;
    
    model DrumBoiler 
      "Complete drum boiler model, including evaporator and supplementary components" 
      
      import Modelica.SIunits.Conversions.*;
      
      Evaporator evaporator annotation (extent=[-50, -20; -30, 0]);
      annotation (
        uses(Modelica_Fluid(version="0.72")),
        Diagram,
        Coordsys(
          extent=[-100, -100; 100, 100],
          grid=[1, 1],
          component=[20, 20]),
        Icon(
          Rectangle(extent=[-100, 100; 100, -100], style(fillColor=7)),
          Text(
            extent=[-151, 165; 138, 102],
            style(fillColor=7, fillPattern=1),
            string="%name"),
          Text(
            extent=[-79, 67; 67, 21],
            string="drum",
            style(
              color=0,
              fillColor=7,
              fillPattern=1)),
          Text(
            extent=[-90, -14; 88, -64],
            string="boiler",
            style(
              color=0,
              fillColor=7,
              fillPattern=1))));
      Modelica.Thermal.HeatTransfer.PrescribedHeatFlow furnace 
        annotation (extent=[-50, -50; -30, -30], rotation=90);
      Modelica.Blocks.Interfaces.RealInput q_F(redeclare type SignalType = Real (
             unit="MW")) 
        annotation (extent=[-109,-55; -100,-65]);
      Modelica.Blocks.Interfaces.RealInput Y_Valve(redeclare type SignalType = 
            Real (unit="1")) 
        annotation (extent=[-109,-85; -100,-75]);
      Valve valve(k = 1.5e-5, redeclare package Medium = 
            Modelica_Media.Water.WaterIF97_ph) 
                    annotation (extent=[44, -20; 64, 0]);
      Modelica_Fluid.Sources.FixedAmbient sink(redeclare package Medium = 
            Modelica_Media.Water.WaterIF97OnePhase_pT) 
        annotation (extent=[80, -20; 100, 0], rotation=180);
      Modelica_Fluid.Sensors.MassFlowRate massFlowRate(redeclare package Medium
          =        Modelica_Media.Water.WaterIF97_ph) 
        annotation (extent=[10, -20; 30, 0], rotation=180);
      Modelica_Fluid.Sensors.Temperature temperature(
          redeclare package Medium = Modelica_Media.Water.WaterIF97_ph) 
        annotation (extent=[10, 60; 30, 80]);
      Modelica_Fluid.Sensors.Pressure pressure(redeclare package Medium = 
                   Modelica_Media.Water.WaterIF97_ph) 
        annotation (extent=[10,24; 30,44]);
      Modelica.Blocks.Continuous.PI controller(T=120, k=10) 
        annotation (extent=[-51,33; -65,47]);
      Modelica_Fluid.Sources.PrescribedMassFlowRate_hX pump(redeclare package 
          Medium = 
            Modelica_Media.Water.WaterIF97_ph, h_ambient=5e5) 
        annotation (extent=[-80, -20; -60, 0]);
      Modelica.Blocks.Math.Feedback feedback 
        annotation (extent=[-26, 30; -46, 50]);
      Modelica.Blocks.Sources.Constant levelSetPoint(k=67) 
        annotation (extent=[-43,60; -30,73]);
      Modelica.Blocks.Interfaces.RealOutput T_S(redeclare type SignalType = 
            Real (unit="degC")) 
        annotation (extent=[100,66; 108,74]);
      Modelica.Blocks.Interfaces.RealOutput p_S(redeclare type SignalType = 
            Real (unit="bar")) 
        annotation (extent=[100,30; 108,38]);
      Modelica.Blocks.Interfaces.RealOutput qm_S(redeclare type SignalType = 
            Modelica.SIunits.MassFlowRate) 
        annotation (extent=[100,6; 108,14],    rotation=0);
      Modelica.Blocks.Interfaces.RealOutput sigma_D 
        annotation (extent=[-24,1; -16,9]);
      Modelica.Blocks.Interfaces.RealOutput V_l(redeclare type SignalType = 
            Modelica.SIunits.Volume) 
        annotation (extent=[-24,16; -16,24]);
      Modelica.Blocks.Math.Gain MW2W(k=1e6) 
        annotation (extent=[-95,-65.5; -85,-54.5]);
      Modelica.Blocks.Math.Gain Pa2bar(k=1e-5) annotation (extent=[37,29; 47,39]);
      Modelica.Thermal.HeatTransfer.Celsius.FromKelvin K2degC 
        annotation (extent=[38,65; 48,75]);
      Modelica.Blocks.Nonlinear.Limiter limiter(uMin=0, uMax=500)
        annotation (extent=[-85,33; -71,47], rotation=180);
    equation 
      connect(furnace.port, evaporator.heatPort) 
        annotation (points=[-40, -30; -40, -21], style(color=42));
      connect(Y_Valve, valve.Y) 
        annotation (points=[-104.5,-80; 54,-80; 54,-17],  style(rgbcolor={0,0,127}));
      connect(evaporator.port_b, temperature.port) annotation (points=[-29,-10;
            -2,-10; -2,51; 20,51; 20,59],               style(color=69));
      connect(evaporator.port_b, pressure.port) annotation (points=[-29,-10; -2,
            -10; -2,15; 20,15; 20,23],              style(color=69));
      connect(evaporator.port_b, massFlowRate.port_a) 
        annotation (points=[-29,-10; -5,-10; -5,-10; 31,-10],
                                               style(color=69));
      connect(massFlowRate.port_b, valve.port_a) 
        annotation (points=[9,-10; 32,-10; 32,-10; 43,-10],
                                               style(color=69));
      connect(valve.port_b, sink.port) annotation (points=[65,-10; 72,-10; 72,
            -10; 79,-10],      style(color=69));
      connect(pump.port, evaporator.port_a) 
        annotation (points=[-59, -10; -51, -10], style(color=69));
      connect(controller.u,feedback.y) 
        annotation (points=[-49.6,40; -45,40], style(rgbcolor={0,0,127}));
      connect(feedback.u2,      evaporator.V) 
        annotation (points=[-36, 32; -36, 1], style(rgbcolor={0,0,127}));
      connect(levelSetPoint.y,feedback.u1)             annotation (points=[
            -29.35,66.5; -22,66.5; -22,40; -28,40],
                                             style(rgbcolor={0,0,127}));
      connect(massFlowRate.m_flow, qm_S) 
        annotation (points=[20,1; 20,10; 104,10],       style(rgbcolor={0,0,127}));
      connect(evaporator.sigma_D, sigma_D) annotation (points=[-29,-5; -26,-5;
            -26,5; -20,5],         style(rgbcolor={0,0,127}));
      connect(evaporator.V, V_l) 
        annotation (points=[-36,1; -36,20; -20,20],    style(rgbcolor={0,0,127}));
      connect(MW2W.y,furnace.Q_flow)       annotation (points=[-84.5,-60; -40,
            -60; -40,-50],     style(rgbcolor={0,0,127}));
      connect(pressure.p, Pa2bar.u) 
        annotation (points=[31,34; 36,34], style(color=74, rgbcolor={0,0,127}));
      connect(Pa2bar.y, p_S) 
        annotation (points=[47.5,34; 104,34],style(color=74, rgbcolor={0,0,127}));
      connect(q_F, MW2W.u) annotation (points=[-104.5,-60; -96,-60],
                                                                   style(color=74,
            rgbcolor={0,0,127}));
      connect(temperature.T, K2degC.Kelvin) annotation (points=[31,70; 37,70],
          style(color=74, rgbcolor={0,0,127}));
      connect(K2degC.Celsius, T_S) annotation (points=[48.5,70; 104,70],style(
            color=74, rgbcolor={0,0,127}));
      connect(controller.y, limiter.u) annotation (points=[-65.7,40; -69.6,40], 
          style(color=74, rgbcolor={0,0,127}));
      connect(limiter.y, pump.m_flow_ambient) annotation (points=[-85.7,40; -90,
            40; -90,-10; -82,-10], style(color=74, rgbcolor={0,0,127}));
    end DrumBoiler;
    
    package WaterPhaseBoundaryIF97 
      "Physical properties for water at phase boundary at boiling and dew curves" 
      
      extends Modelica_Media.Interfaces.PartialMedium(
        mediumName="WaterIF97",
        substanceNames=fill("", 0),
        incompressible=false,
        reducedX=true,
        MassFlowRate(quantity="MassFlowRate.WaterIF97"));
      
      redeclare model extends BaseProperties 
        
      annotation(structurallyIncomplete);
        parameter Integer region=0 "specify region 1 (liquid) or 2 (vapour)";
      equation 
        
        assert(region == 1 or region == 2,
          "WaterPhaseBoundaryIF97 medium model only valid for regions 1 and 2");
        T = Modelica_Media.Water.IF97_Utilities.BaseIF97.Basic.tsat(p);
        if region == 1 then
          d = Modelica_Media.Water.IF97_Utilities.BaseIF97.Regions.rhol_p(p);
          h = Modelica_Media.Water.IF97_Utilities.BaseIF97.Regions.hl_p(p);
        else
          d = Modelica_Media.Water.IF97_Utilities.BaseIF97.Regions.rhov_p(p);
          h = Modelica_Media.Water.IF97_Utilities.BaseIF97.Regions.hv_p(p);
        end if;
        u = h - p/d;
        R = Modelica.Constants.R / Modelica_Media.Water.IF97_Utilities.BaseIF97.data.MH2O;
      end BaseProperties;
    end WaterPhaseBoundaryIF97;
  end Components;
end DrumBoiler;
