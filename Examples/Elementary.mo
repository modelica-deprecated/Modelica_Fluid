  package Elementary 
    "Elementary examples to demonstrate various features of the fluid library"
     
    
    extends Modelica.Icons.Library;
    model SimpleMixing 
      "This example shows the difference of a JunctionVolume and of a FixedComponentVolume"
       
      
      Interfaces.JunctionVolume junctionVolume(
        V=1.e-4, 
        T_start=from_degC(50.0), 
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater)
        annotation (extent=[-10, 40; 10, 20], rotation=180);
      
      import Modelica.SIunits.Conversions.*;
      
      extends Modelica.Icons.Example;
      Sources.FixedAmbient fixedAmbient1(
        p_ambient=from_bar(1.0), 
        T_ambient=from_degC(0), 
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater)
        annotation (extent=[-80, 20; -60, 40]);
      Components.ShortPipe shortPipe1(
        m_dot_nominal=10, 
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater, 
        frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar)
        annotation (extent=[-40, 20; -20, 40]);
      annotation (Diagram, Documentation(info="<html>
<p>
This example shows two ways of merging flows. In the upper part
a volume in the connection port is modeled. At the initial time the
fluid in this volume has a higher temperature as the fluid
flowing form the ambient. Therefore, im shortPipe the temperature
will need some time to arrive at the mixing temperature.
</p>
<p>
In the lower part the volume in the connection port is neglected
and ideal mixing takes place. Since there is a steady state mass
flow and no (existing) delay mechanism is modeled, the mixing
is instantaneous. The temperature is the same as in the upper part.
</p>
</html>
"));
      Components.ShortPipe shortPipe3(
        m_dot_nominal=10, 
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater, 
        frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar)
        annotation (extent=[20, 20; 40, 40]);
      Sources.FixedAmbient fixedAmbient3(
        T_ambient=from_degC(20), 
        p_ambient=from_bar(0.95), 
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater)
        annotation (extent=[80, 20; 60, 40]);
      Sources.FixedAmbient fixedAmbient2(
        p_ambient=from_bar(1.0), 
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater, 
        T_ambient=from_degC(10)) annotation (extent=[-80, 60; -60, 80]);
      Components.ShortPipe shortPipe2(
        m_dot_nominal=10, 
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater, 
        frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar)
        annotation (extent=[-40, 60; -20, 80]);
      Sources.FixedAmbient fixedAmbient4(
        p_ambient=from_bar(1.0), 
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater, 
        T_ambient=from_degC(1)) annotation (extent=[-80, -80; -60, -60]);
      Components.ShortPipe shortPipe4(
        m_dot_nominal=10, 
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater, 
        frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar)
        annotation (extent=[-40, -80; -20, -60]);
      Components.ShortPipe shortPipe6(
        m_dot_nominal=10, 
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater, 
        frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar)
        annotation (extent=[20, -80; 40, -60]);
      Sources.FixedAmbient fixedAmbient6(
        T_ambient=from_degC(20), 
        p_ambient=from_bar(0.95), 
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater)
        annotation (extent=[80, -80; 60, -60]);
      Sources.FixedAmbient fixedAmbient5(
        p_ambient=from_bar(1.0), 
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater, 
        T_ambient=from_degC(10)) annotation (extent=[-80, -40; -60, -20]);
      Components.ShortPipe shortPipe5(
        m_dot_nominal=10, 
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater, 
        frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar)
        annotation (extent=[-40, -40; -20, -20]);
    equation 
      connect(fixedAmbient1.port, shortPipe1.port_a)
        annotation (points=[-59, 30; -41, 30], style(color=69));
      connect(shortPipe3.port_b, fixedAmbient3.port)
        annotation (points=[41, 30; 59, 30], style(color=69));
      connect(shortPipe3.port_a, junctionVolume.port)
        annotation (points=[19, 30; 0, 30], style(color=69));
      connect(fixedAmbient2.port, shortPipe2.port_a)
        annotation (points=[-59, 70; -41, 70], style(color=69));
      connect(fixedAmbient4.port, shortPipe4.port_a)
        annotation (points=[-59, -70; -41, -70], style(color=69));
      connect(shortPipe6.port_b, fixedAmbient6.port)
        annotation (points=[41, -70; 59, -70], style(color=69));
      connect(fixedAmbient5.port, shortPipe5.port_a)
        annotation (points=[-59, -30; -41, -30], style(color=69));
      connect(shortPipe5.port_b, shortPipe6.port_a) annotation (points=[-19, -30; 
             0, -30; 0, -70; 19, -70], style(color=69));
      connect(shortPipe4.port_b, shortPipe6.port_a)
        annotation (points=[-19, -70; 19, -70], style(color=69));
      connect(shortPipe2.port_b, junctionVolume.port)
        annotation (points=[-19, 70; 0, 70; 0, 30], style(color=69));
      connect(shortPipe1.port_b, junctionVolume.port)
        annotation (points=[-19, 30; 0, 30], style(color=69));
    end SimpleMixing;
    
  end Elementary;
