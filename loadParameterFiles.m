function params = loadParameterFiles(ID)
    switch ID
        case 1 
            params = loadParameters_YBCO610_Bonding2_U0_0_0;
        case 2 
            params = loadParameters_YBCO655_Bonding2_U0_0_0;
        case 3 
            params = loadParameters_YBCO679_Bonding2_U0_0_0;
        case 4
            params = loadParameters_YBCO610_Bonding2_U0_1_0;
        case 5 
            params = loadParameters_YBCO655_Bonding2_U0_1_0;
        case 6 
            params = loadParameters_YBCO679_Bonding2_U0_1_0;
        case 7 
            params = loadParameters_YBCO699_Bonding2_U0_1_0;
        case 8 
            params = loadParameters_YBCOCa_Bonding2_U0_1_0;
        case 9 
            params = loadParameters_YBCO610_Bonding2_U0_0_5;
        case 10 
            params = loadParameters_YBCO655_Bonding2_U0_0_5;
        case 11
            params = loadParameters_YBCO679_Bonding2_U0_0_5;
        case 12
            params = loadParameters_YBCO699_Bonding2_U0_0_5;
        case 13 
            params = loadParameters_YBCOCa_Bonding2_U0_0_5;
        case 14
            params = loadParameters_YBCO655_Bonding2_U0_1_0__All_momenta;
        case 15
            params = loadParameters_YBCO610_Bonding2_U0_1_5;
        case 16
            params = loadParameters_YBCO655_Bonding2_U0_1_5;
        case 17
            params = loadParameters_YBCO679_Bonding2_U0_1_5;
        case 18
            params = loadParameters_YBCO699_Bonding2_U0_1_5;
        case 19
            params = loadParameters_YBCOCa_Bonding2_U0_1_5;
        case 20
            params = loadParameters_YBCO610_Bonding2_U0_2_0;
        case 21
            params = loadParameters_YBCO655_Bonding2_U0_2_0;
        case 22
            params = loadParameters_YBCO679_Bonding2_U0_2_0;
        case 23
            params = loadParameters_YBCO699_Bonding2_U0_2_0;
        case 24
            params = loadParameters_YBCOCa_Bonding2_U0_2_0;
	case 25
	    params = loadParameters_YBCO635_Bonding2_U0_1_0;
    end
end
