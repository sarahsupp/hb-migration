CREATE TABLE access_key_lut (
	ak_id serial NOT NULL PRIMARY KEY,
	access_key bigint NOT NULL UNIQUE,
	date_time timestamp,
	tracks text,
	xml text,
	status text,
	ann_name text,
	spp text
	);