CREATE TABLE access_key_lut (
	ak_id serial NOT NULL PRIMARY KEY,
	access_key numeric NOT NULL,
	date_time text,
	tracks text,
	xml text,
	status text
	);