{
	/*
	#include <thread>
	#include <chrono>
	std::vector<std::pair<std::string,std::string>> stuff{{"raw","275_A"},
	                                                     {"simple","275_A"},
	                                                     {"complex","275_A"},
	                                                     {"raw","275_B"},
	                                                     {"simple","275_B"},
	                                                     {"complex","275_B"}};
	TCanvas c1;
	for(auto&& things : stuff){
		std::string name=things.first;
		std::string calib_ver=things.second;
	*/
		std::string name="275_A";
		std::string calib_ver="complex";
		TCanvas c1;
		TF1 calib_curve("calib", "pol6", 0, 0.25);
		//calib_curve.SetNpx(1000);
		if(calib_ver=="raw" && name=="275_A"){
			//std::vector<double> calib_coefficients{0.00517657, 2.86027, -7.36815, -19.955, 78.0564, 491.547, -1766.58}; // jul08
			std::vector<double> calib_coefficients{0.00141186, 2.85759, -7.43248, -18.2004, 78.8485, 472.124, -1814.29};
			calib_curve.SetParameters(calib_coefficients.data());
		} else if(calib_ver=="simple" && name=="275_A"){
			//std::vector<double> calib_coefficients{0.00100964, 3.03952, -19.4798, 77.005, 27.6594, -1608.27, 4177.23}; // jul08
			std::vector<double> calib_coefficients{0.000901656, 2.7328, -11.972, 19.2158, 70.2304, -429.644, 743.848}; // jul12
			calib_curve.SetParameters(calib_coefficients.data());
		} else if(calib_ver=="complex" && name=="275_A"){
			//std::vector<double> calib_coefficients{0.00532017, 2.62183, -8.00979, -18.5983, 91.0425, 516.762, -1989.37}; // jul08
			std::vector<double> calib_coefficients{0.00256741, 2.53608, -5.59115, -34.2931, 86.1831, 887.119, -2934}; // jul12
			calib_curve.SetParameters(calib_coefficients.data());
		} else if(calib_ver=="raw" && name=="275_B"){
			//std::vector<double> calib_coefficients{0.0121442, 2.90637, -7.30414, -19.9099, 77.1875, 485.482, -1764.29}; // jul08
			std::vector<double> calib_coefficients{0.00151493, 2.67998, -6.92564, -16.6946, 75.3407, 444.124, -1731.61}; // jul12
			calib_curve.SetParameters(calib_coefficients.data());
		} else if(calib_ver=="simple" && name=="275_B"){
			//std::vector<double> calib_coefficients{0.00958699, 3.01173, -17.3047, 56.2216, 55.7341, -1232.64, 2979.2}; // jul08
			std::vector<double> calib_coefficients{0.000557654, 2.58452, -12.3119, 31.4942, 55.7869, -735.182, 1678.52}; // jul12
			calib_curve.SetParameters(calib_coefficients.data());
		} else if(calib_ver=="complex" && name=="275_B"){
			//std::vector<double> calib_coefficients{0.00603398, 2.88702, -9.77863, -18.9483, 116.117, 594.401, -2502.56}; // jul08
			std::vector<double> calib_coefficients{0.00197077, 2.35596, -5.21091, -30.7562, 100.809, 708.71, -2614.03};  // jul12
			calib_curve.SetParameters(calib_coefficients.data());
		}
		calib_curve.Draw();
		/*
		//if(things==stuff.front()) calib_curve.DrawClone();
		//else calib_curve.DrawClone("same");
		std::cout<<"curve: "<<name<<", "<<calib_ver<<std::endl;
		calib_curve.Draw();
		c1.Modified();
		c1.Update();
		gSystem->ProcessEvents();
		std::this_thread::sleep_for(std::chrono::milliseconds(3000));
	}*/
}
