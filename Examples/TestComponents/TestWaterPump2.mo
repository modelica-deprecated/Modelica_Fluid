model TestWaterPump2 
  import PC = Modelica_Fluid.Components.PumpCharacteristics;
  function pumpConsChar2 = PC.QuadraticPowerCharacteristic (
    q_nom={0,0.001,0.0015}, P_nom={550,650,800}, n0=1200);
  extends TestWaterPump(Pump1(M=0.1,redeclare function powerCharacteristic = 
          pumpConsChar2));
  
  annotation (Diagram);
end TestWaterPump2;
