within FluidSandbox;
package Turbomachinery "Compressor, turbine models and the like" 
  extends Icons.VariantLibrary;
  model TurbineStage "Turbine stage according to Stodola's law" 
    
    extends Interfaces.PartialComponent;
    extends FluidInterface.PartialTransportIsentropic;
    
    parameter Real K_t=0.001 "Stodola's machine constant";
    
    Modelica.Mechanics.Rotational.Interfaces.Flange_a flange 
    annotation (extent=[90,20; 110,40]);
  equation 
    port_a.m_flow = K_t*Modelica_Fluid.Utilities.regRoot((port_a.p^2 - port_b.p^2)/
      medium_designDirection.T) "Stodola's law";
    0 = P_mechanical + flange.tau*der(flange.phi);
    
    annotation (Icon(Polygon(points=[-100,100; 100,150; 100,-50; -100,0; -100,
              100],
            style(
            color=3,
            rgbcolor={0,0,255},
            fillColor=46,
            rgbfillColor={216,62,1}))),      Diagram);
  end TurbineStage;
end Turbomachinery;
