within Modelica_Fluid.Test.TestComponents.Pumps;
model TestWaterPump2
  extends TestWaterPumpDefault(Pump1(M=0.1,
      usePowerCharacteristic=true,
      redeclare function powerCharacteristic = 
          Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics.quadraticPower (
            q_nom={0,0.001,0.0015}, W_nom={550,650,800})),
    Source(T=ambient.default_T_ambient),
    Sink(T=ambient.default_T_ambient));

  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                      graphics),
    experiment(StopTime=10),
    experimentSetupOutput);

end TestWaterPump2;
