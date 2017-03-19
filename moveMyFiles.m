function moveMyFiles()
	delete('*.err');
	delete('*.out');
	movefile('*IRIXS*.mat', '../IRIXS_DATA_15_24/');
	movefile('*XAS*.mat', '../XAS_DATA_15_24/');
end
