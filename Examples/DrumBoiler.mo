package DrumBoiler 
  "Drum boiler example, see Franke, Rode, Krueger: On-line Optimization of Drum Boiler Startup, 3rd International Modelica Conference, Linkoping, 2003" 
  
  model DrumBoilerSimulation "Simulate start-up of DrumBoiler" 
    extends Modelica.Icons.Example;
    Components.DrumBoiler drumBoiler annotation (extent=[-20, -40; 40, 20]);
    Modelica.Blocks.Sources.TimeTable q_F_Tab(table=[0, 0; 3600, 400; 7210,
          400]) annotation (extent=[-80, 2; -60, 22]);
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
    connect(q_F_Tab.y,       drumBoiler.q_F) annotation (points=[-59,12; -40,12; 
          -40,-16; -25.7,-16],        style(color=3));
    connect(Y_Valve_Tab.y,       drumBoiler.Y_Valve) annotation (points=[-59,
          -30; -44,-30; -44,-34; -25.7,-34],     style(
        color=3,
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
      sigma_D           = 60*der(T_D);
      
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
    
    model MassFlowSource_h 
      "Ideal pump that produces a mass flow rate from a large reservoir defined by input signal" 
      
      import Modelica_Fluid.*;
      import Modelica_Fluid.Interfaces.*;
      import Modelica_Media.Interfaces.*;
      
      import SI = Modelica.SIunits;
      import Modelica.SIunits.Conversions.*;
      
      FluidPort_b port(redeclare model Medium = Medium) 
        annotation (extent=[100, -10; 120, 10], rotation=0);
      replaceable package Medium = PartialMedium extends PartialMedium 
        "Medium in the component"   annotation (choicesAllMatching=true);
      
      parameter SI.SpecificEnthalpy h_ambient=5e5 "Ambient enthalphy";
      parameter SI.Pressure p_start=from_bar(1.0) 
        "|Initialization|| Initial pressure";
      annotation (
        Coordsys(
          extent=[-100, -100; 100, 100],
          grid=[2, 2],
          component=[20, 20]),
        Icon(
          Rectangle(extent=[20, 60; 100, -60], style(
              color=0,
              gradient=2,
              fillColor=8)),
          Rectangle(extent=[38, 40; 100, -40], style(
              color=69,
              gradient=2,
              fillColor=69)),
          Ellipse(extent=[-100, 80; 60, -80], style(fillColor=7)),
          Polygon(points=[-60, 70; 60, 0; -60, -68; -60, 70], style(color=73,
                 fillColor=73)),
          Text(
            extent=[-54, 32; 16, -30],
            string="m'",
            style(color=41, fillColor=41)),
          Text(extent=[-142, 142; 156, 88], string="%name"),
          Text(
            extent=[-126, -86; 146, -116],
            style(color=0),
            string="%h_ambient")),
        Window(
          x=0.45,
          y=0.01,
          width=0.44,
          height=0.65),
        Diagram);
      Modelica.Blocks.Interfaces.RealInput m_flow(
                                              redeclare type SignalType = 
            SI.MassFlowRate) 
        "Mass flow rate from an infinite reservoir in to the port as signal" 
        annotation (extent=[-140, -20; -100, 20]);
    equation 
      port.m_flow = -noEvent(max(m_flow,           0));
      port.H_flow = semiLinear(port.m_flow, port.h, h_ambient);
    end MassFlowSource_h;
    
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
      Modelica.Blocks.Interfaces.RealInput q_F 
        annotation (extent=[-139, 0; -99, -40]);
      Modelica.Blocks.Interfaces.RealInput Y_Valve 
        annotation (extent=[-139, -100; -99, -60]);
      Valve valve(                                                           k=
            1.5e-5, redeclare package Medium = 
            Modelica_Media.Water.StandardWater) 
                    annotation (extent=[44, -20; 64, 0]);
      Modelica_Fluid.Sources.FixedAmbient sink(redeclare package Medium = 
            Modelica_Media.Water.OnePhaseWater) 
        annotation (extent=[80, -20; 100, 0], rotation=180);
      Modelica_Fluid.Sensors.MassFlowRate massFlowRate(redeclare package Medium
          =        Modelica_Media.Water.WaterIF97) 
        annotation (extent=[10, -20; 30, 0]);
      Modelica_Fluid.Sensors.Temperature temperature(      signalUnit="degC",
          redeclare package Medium = Modelica_Media.Water.StandardWater) 
        annotation (extent=[10, 60; 30, 80]);
      Modelica_Fluid.Sensors.Pressure pressure(redeclare package Medium = 
                   Modelica_Media.Water.WaterIF97, signalUnit="bar") 
        annotation (extent=[10, 20; 30, 40]);
      Modelica.Blocks.Continuous.PI controller(T=120, k=10) 
        annotation (extent=[-60, 30; -80, 50]);
      MassFlowSource_h pump(redeclare package Medium = 
            Modelica_Media.Water.WaterIF97) 
        annotation (extent=[-80, -20; -60, 0]);
      Modelica.Blocks.Math.Feedback feedback 
        annotation (extent=[-26, 30; -46, 50]);
      Modelica.Blocks.Sources.Constant levelSetPoint(k=67) 
        annotation (extent=[-46, 60; -26, 80]);
    protected 
      Modelica.Blocks.Interfaces.RealOutput T_S 
        annotation (extent=[34, 66; 42, 74]);
      Modelica.Blocks.Interfaces.RealOutput p_S 
        annotation (extent=[34, 26; 42, 34]);
      Modelica.Blocks.Interfaces.RealOutput qm_S 
        annotation (extent=[34, -34; 42, -26], rotation=0);
      Modelica.Blocks.Interfaces.RealOutput sigma_D 
        annotation (extent=[-24, 6; -16, 14]);
      Modelica.Blocks.Interfaces.RealOutput V_l 
        annotation (extent=[-24, 24; -16, 32]);
    public 
      Modelica.Blocks.Math.Gain MW2W(k=1e6) 
        annotation (extent=[-70, -69; -50, -50]);
    equation 
      connect(furnace.port, evaporator.heatPort) 
        annotation (points=[-40, -30; -40, -21], style(color=42));
      connect(Y_Valve, valve.Y) 
        annotation (points=[-119, -80; 54, -80; 54, -17], style(color=3));
      connect(evaporator.port_b, temperature.port) annotation (points=[
            -29, -10; -2, -10; -2, 50; 20, 50; 20, 59], style(color=69));
      connect(evaporator.port_b, pressure.port) annotation (points=[-29,
             -10; -2, -10; -2, 10; 20, 10; 20, 19], style(color=69));
      connect(evaporator.port_b, massFlowRate.port_a) 
        annotation (points=[-29, -10; 9, -10], style(color=69));
      connect(massFlowRate.port_b, valve.port_a) 
        annotation (points=[31, -10; 43, -10], style(color=69));
      connect(valve.port_b, sink.port) annotation (points=[65,-10; 72,-10; 72,
            -10; 79,-10],      style(color=69));
      connect(pump.port, evaporator.port_a) 
        annotation (points=[-59, -10; -51, -10], style(color=69));
      connect(controller.u,feedback.y) 
        annotation (points=[-58, 40; -45, 40], style(color=3));
      connect(feedback.u2,      evaporator.V) 
        annotation (points=[-36, 32; -36, 1], style(color=3));
      connect(levelSetPoint.y,feedback.u1)             annotation (points=[-25,
             70; -20, 70; -20, 40; -28, 40], style(color=3));
      connect(pressure.p, p_S) 
        annotation (points=[31,30; 38,30],   style(color=3));
      connect(temperature.T, T_S) 
        annotation (points=[31,70; 38,70],   style(color=3));
      connect(massFlowRate.m_flow, qm_S) 
        annotation (points=[20,-21; 20,-30; 38,-30],    style(color=3));
      connect(evaporator.sigma_D, sigma_D) annotation (points=[-29, -5; -26,
            -5; -26, 10; -20, 10], style(color=3));
      connect(evaporator.V, V_l) 
        annotation (points=[-36, 1; -36, 28; -20, 28], style(color=3));
      connect(controller.y,       pump.m_flow) annotation (points=[-81, 40; -90,
             40; -90, -10; -82, -10], style(color=3));
      connect(q_F,MW2W.u)       annotation (points=[-119, -20; -90, -20; -90,
             -59.5; -72, -59.5], style(color=3));
      connect(MW2W.y,furnace.Q_flow)       annotation (points=[-49, -59.5; -40,
             -59.5; -40, -50], style(color=3));
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
        T = Modelica_Media.Water.IF97.BaseIF97.Basic.tsat(p);
        if region == 1 then
          d = Modelica_Media.Water.IF97.BaseIF97.Regions.rhol_p(p);
          h = Modelica_Media.Water.IF97.BaseIF97.Regions.hl_p(p);
        else
          d = Modelica_Media.Water.IF97.BaseIF97.Regions.rhov_p(p);
          h = Modelica_Media.Water.IF97.BaseIF97.Regions.hv_p(p);
        end if;
        u = h - p/d;
        R = 287.0; // data.R // Modelica.Constants.R/data.MM;
      end BaseProperties;
    end WaterPhaseBoundaryIF97;
  end Components;
end DrumBoiler;
