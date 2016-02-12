#!/usr/bin/perl
use strict;
use FindBin qw ($Bin);
use Getopt::Long;
use SVG;

my %opt;
GetOptions(\%opt,"table:s","help:s");

my $help=<<USAGE;
perl $0 --table ../input/QTL.regions.list
USAGE

if(defined $opt{help} or keys %opt < 1){
        die  $help ;
}



my $svg=SVG->new(width=>1400,height=>1000);
my $startw=50; 
my $endw  =1300;
my $starth=100;
my $endh  =500;
my $scale  = 500/4500;

print "Plot type1\n";
plot_deletion_two_mping($svg, $startw, $starth, 'Type 1: Deletions between two nearby mPing insertions', 'RIL_SV_type1.txt', $scale);
print "Plot type2\n";
plot_deletion_two_mping($svg, $startw, $starth+500, 'Type 2: Deletions from one mPing insertions to flanking sequence', 'RIL_SV_type2.txt', $scale);
print "Plot type3\n";
plot_deletion_two_mping($svg, $startw+600, $starth, 'Type 3: Deletions accompanied by rearrangement of flanking sequence', 'RIL_SV_type3.txt', $scale);
#plot_deletion_one_mping($svg, $startw+500, $starth, 'Type 2: Deletions from one mPing insertions to flanking sequence');
#plot_deletion_one_mping($svg, $startw, $starth+350, 'Type 3: Deletions accompanied by rearrangement of flanking sequence');
#plot_deletion_one_mping()
#plot_deletion_rearrangement()

plot_scale_bar($svg, $scale);

writesvg("mPing_SV_diagram.svg", $svg);

###sub functions
sub plot_deletion_two_mping
{
my ($svg, $startx, $starty, $title, $plot_inf, $scale) = @_;
my $plot_data = read_inf($plot_inf);
plot_text($svg, $startx+15, $starty-50, 'start', 15, $title);

foreach my $rank (sort {$a <=> $b} keys %$plot_data){
    print "SV$rank\n";
    my $x = $startx;
    my $y = $starty + ($rank-1)*150;
    plot_diagram_simple($svg, $plot_data->{$rank}, $x, $y, $scale);
}
}

##############################sub plot scale bar and mping box
sub plot_scale_bar
{
my ($svg, $scale) = @_;

#scale bar
my $x1 = 800;
my $x2 = $x1 + $scale*1000;
my $y1 = 600;
my $y2 = $y1;
my $line=$svg->line(
           x1=>$x1, y1=>$y1,
           x2=>$x2, y2=>$y2,
           style=>{
               stroke=>'black'
           }
       );
my $scale_text =$svg->text(
           x=>$x1 + 120, y=>$y1,
           style=>{
               fontsize=>'4','text-anchor'=>'start','stroke-width'=>1
           }
       )->cdata('1 kb');

#mping on direct diretion
my $mping_plus_y = 630;
my $mping_plus_x1 = 800;
my $mping_plus_x2 = $mping_plus_x1 + $scale*500;
mping_box($svg, $mping_plus_x1+50, $mping_plus_x2+50, "+", $mping_plus_y);
my $mping_plus =$svg->text(
           x=>$x1 + 120, y=>$mping_plus_y,
           style=>{
               'font-size'=>'15','text-anchor'=>'start', 'font-style'=>'italic', 'stroke-width'=>1
           }
       )->cdata('mPing on plus strand');

#mping in reverse direction
my $mping_minus_y = 680;
mping_box($svg, $mping_plus_x1+50, $mping_plus_x2+50, "-", $mping_minus_y);
my $mping_minus =$svg->text(
          x=>$x1 + 120, y=>$mping_minus_y,
          style=>{
               'font-size'=>'15','text-anchor'=>'start', 'font-style'=>'italic', 'stroke-width'=>1
          }
       )->cdata('mPing on minus strand');


}

################################### sub for plot diagram of one SV#####################
sub plot_diagram_simple
{
my ($svg, $data, $x, $y, $scale) = @_;

#header for this plot
plot_text($svg, $x-25, $y-25, 'start', '20', $data->{'header'});

#upper plot: HEG4
plot_text($svg, $x, $y, 'start', '15', $data->{'up_header'});
#plot_chromosome_line()
#plot_zoom_in_mping($svg, $x, $y, 'up', $data) if exists $data->{'up_zoom_in_mping'};
#plot_zoom_in($svg, $x, $y, 'up', $data) if exists $data->{'up_zoom_in'};
my $up_fragment_hash = plot_fragment_arrow($svg, $x+50, $y-10, $data->{'up_segment'}, $data->{'up_cut'}, $data->{'up_start'}, $scale, -10);
my $up_mping_hash    = plot_mping($svg, $x+50, $y-10, $data->{'up_zoom_in_mping'}, $data->{'up_cut'}, $data->{'up_start'}, $scale);

#lower plot: RILs
plot_text($svg, $x, $y+70, 'start', '15', $data->{'down_header'});
#plot_chromosome_line() 
#plot_zoom_in_mping($svg, $x, $y+70, 'down', $data) if exists $data->{'down_zoom_in_mping'};
#plot_zoom_in($svg, $x, $y+70, 'down', $data) if exists $data->{'down_zoom_in'}; 
my $down_fragment_hash = plot_fragment_arrow($svg, $x+50, $y+70, $data->{'down_segment'}, $data->{'down_cut'}, $data->{'down_start'}, $scale, 15);
my $down_mping_hash    = plot_mping($svg, $x+50, $y+70, $data->{'down_zoom_in_mping'}, $data->{'down_cut'}, $data->{'down_start'}, $scale);

#comparison plot
#my $up_hash = merge_hash();
@$up_fragment_hash{keys %$up_mping_hash} = values %$up_mping_hash;
@$down_fragment_hash{keys %$down_mping_hash} = values %$down_mping_hash;
print "fragment position in figure\n";
print "up\n";
foreach my $t (keys %$up_fragment_hash){
   print "$t\t$up_fragment_hash->{$t}->{'start'}\t$up_fragment_hash->{$t}->{'end'}\n";
}
print "down\n";
foreach my $t (keys %$down_fragment_hash){
   print "$t\t$down_fragment_hash->{$t}->{'start'}\t$down_fragment_hash->{$t}->{'end'}\n";
}
#plot_matches($svg, $x + 50, $y, $x + 50, $y+60, $data);
plot_matches_by_labels($svg, $x + 50, $y, $x + 50, $y+60, $up_fragment_hash, $down_fragment_hash, $data);
}


#plot fragment with arrow on chromosome
#up_segment=a:28723476:28724476:arrow-right;bc:28724476:28732623:arrow-both;d:28732623-28733623:arrow-left
sub plot_fragment_arrow
{
my ($svg, $x, $y, $segment, $cut, $start, $scale, $label_y_adj) = @_;
my %frag_hash;
my @fragment = split(';', $segment);
for (my $i=0; $i<@fragment; $i++){
    my ($label, $f_start, $f_end, $arrow) = split(':', $fragment[$i]);
    my $f_start_rl = convert_position($cut, $start, $f_start, $scale, $x);
    my $f_end_rl = convert_position($cut, $start, $f_end, $scale, $x);
    $frag_hash{$label}{'start'} = $f_start_rl;
    $frag_hash{$label}{'end'} = $f_end_rl;
    print "$label\t$f_start_rl\t$f_end_rl\t$arrow\n";
    if ($label eq 'hi'){
        $label = 'de';
    }
    #draw arrow on left of segment
    if ($arrow eq 'arrow-left'){
       #line and arrow
       my $line=$svg->line(
           x1=>$f_start_rl+10, y1=>$y,
           x2=>$f_end_rl, y2=>$y,
           style=>{stroke=>'black', 'marker-start'=>'url(#left_arrow)'}
       );
       #label
       my $text=$svg->text(
           x=>$f_start_rl+15, y=>$y+$label_y_adj,
               style=>{
                  fontsize=>'4','text-anchor'=>'start','stroke-width'=>1
               }
       )->cdata($label);
    #draw arrow on right of segment
    }elsif($arrow eq 'arrow-right'){
       #line and arrow
       my $line=$svg->line(
           x1=>$f_start_rl, y1=>$y,
           x2=>$f_end_rl-10, y2=>$y,
           style=>{stroke=>'black', 'marker-end'=>'url(#right_arrow)'}
       );
       #label
       my $text=$svg->text(
           x=>$f_end_rl-25, y=>$y+$label_y_adj,
              style=>{            
                 fontsize=>'4','text-anchor'=>'start','stroke-width'=>1
              }
       )->cdata($label);

    #draw arrow on both side of segment
    }elsif($arrow eq 'arrow-both'){
       #line and arrow
       my $line=$svg->line(
           x1=>$f_start_rl+10, y1=>$y,
           x2=>$f_end_rl-10, y2=>$y,
           style=>{stroke=>'black', 'marker-start'=>'url(#left_arrow)', 'marker-end'=>'url(#right_arrow)'}
       );
       #label on start
       my $text_s=$svg->text(
           x=>$f_start_rl+15, y=>$y+$label_y_adj,
           style=>{              
               'font-size'=>'15','text-anchor'=>'start','stroke-width'=>1
           }
       )->cdata(substr($label, 0, 1));
       #label on end
       my $text_e=$svg->text(
           x=>$f_end_rl-25, y=>$y+$label_y_adj,
           style=>{            
               'font-size'=>'15','text-anchor'=>'start','stroke-width'=>1
           }
       )->cdata(substr($label, 1, 1));
     
    }elsif($arrow eq 'dashed-line'){
       #dashed line, which are case have not figured out
       my $line=$svg->line(
           x1=>$f_start_rl+10, y1=>$y,
           x2=>$f_end_rl-10, y2=>$y,
           style=>{stroke=>'black', 'stroke-dasharray'=>"5,5", 'fill'=>"none"}
       );
    }elsif($arrow eq 'short-arrow'){
       my $line=$svg->line(
           x1=>$f_start_rl+5, y1=>$y,
           x2=>$f_end_rl-5, y2=>$y,
           #style=>{stroke=>'black'}
           style=>{stroke=>'black', 'marker-start'=>'url(#left_arrow)', 'marker-end'=>'url(#right_arrow)'}
       );
       my $text=$svg->text(
           x=>$f_end_rl-10, y=>$y+$label_y_adj,
           style=>{            
               'font-size'=>'15','text-anchor'=>'start','stroke-width'=>1
           }
       )->cdata($label);
    }
    #plot two lines at cut position on chromosome
    plot_cut_lines($svg, $cut, $scale, $start, $x, $y);
}
return \%frag_hash;
}

sub plot_cut_lines
{
my ($svg, $cut, $scale, $start, $x, $y) = @_;
my ($cut_pos, $cut_length) = split(',', $cut);
my $cut_pos_rl = convert_position($cut, $start, $cut_pos, $scale, $x);

my $line1=$svg->line(
       x1=>$cut_pos_rl-2, y1=>$y+4,
       x2=>$cut_pos_rl+2, y2=>$y-4,
       style=>{stroke=>'black'}
    );
my $line2=$svg->line(
       x1=>$cut_pos_rl+3-2, y1=>$y+4,
       x2=>$cut_pos_rl+3+2, y2=>$y-4,
       style=>{stroke=>'black'}
    );

my $length = sprintf('%3d kb', $cut_length/1000);
my $text=$svg->text(
       x=>$cut_pos_rl-15, y=>$y-20,
       style=>{
            'font-size'=>'15','text-anchor'=>'start', 'stroke-width'=>0.1
       }
    )->cdata($length);

}

################################### sub for plot mping ################################
sub plot_mping
{
my ($svg, $x, $y, $mping, $cut, $start, $scale) = @_;
my %mping_hash;
my @mpings = split(';', $mping);
for (my $i=0; $i<@mpings; $i++){
    my ($mping_id, $m_start, $m_strand, $m_left_seq, $m_right_seq) = split(',', $mpings[$i]);
    my $m_end = $m_start + 500;
    my $m_start_rl = convert_position($cut, $start, $m_start, $scale, $x);
    my $m_end_rl = convert_position($cut, $start, $m_end, $scale, $x);
    $mping_hash{$mping_id}{'start'} = $m_start_rl;
    $mping_hash{$mping_id}{'end'} = $m_end_rl;
    mping_box($svg, $m_start_rl, $m_end_rl, $m_strand, $y);
=cut
    my $mping_box = $svg->rectangle(
       x => $m_start_rl, y => $y-10,
       width => $m_end_rl - $m_start_rl + 1, height => 25,
       style=>{
             fill=>'orange',
             'fill-opacity' => 0.7
       }
    );
    if ($m_strand eq '+'){
        print "mping on plus strand\n";
        my $text_plus=$svg->text(
           x=>$m_start_rl, y=>$y+12,
              style=>{            
                 'font-size'=>'15', 'font-style'=>'italic','text-anchor'=>'start','stroke-width'=>1
              }
        )->cdata('mPing');
        my $text=$svg->text(
           x=>$m_start_rl, y=>$y,
              style=>{            
                 'font-size'=>'15', 'font-style'=>'italic','text-anchor'=>'start','stroke-width'=>1
              }
        )->cdata('5\'->3\'');
    }else{
        print "mping on minus strand\n"; 
        my $text_minus=$svg->text(
           x=>$m_end_rl, y=>$y-2,
           style=>{            
                 'font-size'=>'15', 'font-style'=>'italic','text-anchor'=>'start','stroke-width'=>1
           },
           transform => "rotate(180 $m_end_rl, $y)"
        )->cdata('mPing');
        my $text_53=$svg->text(
           x=>$m_start_rl, y=>$y,
              style=>{            
                 'font-size'=>'15', 'font-style'=>'italic','text-anchor'=>'start','stroke-width'=>1
              }
        )->cdata('3\'<-5\'');

    }
=cut
}
return \%mping_hash;
}

sub mping_box
{
my ($svg, $m_start_rl, $m_end_rl, $m_strand, $y) = @_;
    my $mping_box = $svg->rectangle(
       x => $m_start_rl, y => $y-10,
       width => $m_end_rl - $m_start_rl + 1, height => 25,
       style=>{
             fill=>'orange',
             'fill-opacity' => 0.7
       }
    );
    if ($m_strand eq '+'){
        print "mping on plus strand\n";
        my $text_plus=$svg->text(
           x=>$m_start_rl, y=>$y+12,
              style=>{            
                 'font-size'=>'15', 'font-style'=>'italic','text-anchor'=>'start','stroke-width'=>1
              }
        )->cdata('mPing');
        my $text=$svg->text(
           x=>$m_start_rl, y=>$y,
              style=>{            
                 'font-size'=>'15', 'font-style'=>'italic','text-anchor'=>'start','stroke-width'=>1
              }
        )->cdata('5\'->3\'');
    }else{
        print "mping on minus strand\n"; 
        my $text_minus=$svg->text(
           x=>$m_end_rl, y=>$y-2,
           style=>{            
                 'font-size'=>'15', 'font-style'=>'italic','text-anchor'=>'start','stroke-width'=>1
           },
           transform => "rotate(180 $m_end_rl, $y)"
        )->cdata('mPing');
        my $text_53=$svg->text(
           x=>$m_start_rl, y=>$y,
              style=>{            
                 'font-size'=>'15', 'font-style'=>'italic','text-anchor'=>'start','stroke-width'=>1
              }
        )->cdata('3\'<-5\'');

    }
}

################################### sub for plot_strain name ##########################
sub plot_text
{
my ($svg, $x, $y, $anchor, $size, $strain) =@_;
my $strain_name=$svg->text(
               x=>$x, y=>$y,
               style=>{
                    'font-size'=>$size,'text-anchor'=>$anchor,'stroke-width'=>1
               }
    )->cdata($strain);
}

################################## sub for plot match ################################
sub plot_matches_by_labels
{
my ($svg, $x1, $y1, $x2, $y2, $up_hash, $down_hash, $data) = @_;
my @matches = split(';', $data->{'match'});
for (my $i=0; $i<@matches; $i++){
    print "$matches[$i]\n";
    my ($up_label, $down_label, $strand, $color) = split(':', $matches[$i]);
    print "matching pairs: $up_label, $down_label, $strand, $color\n";
    my $match_up_x1 = $up_hash->{$up_label}->{'start'};
    my $match_up_x2 = $up_hash->{$up_label}->{'end'};
    my $match_down_x1 = $down_hash->{$down_label}->{'start'};
    my $match_down_x2 = $down_hash->{$down_label}->{'end'};
    print "match x: $match_up_x1\t$match_up_x2\t$match_down_x1\t$match_down_x2\n";
    print "match y: $y1\t$y1\t$y2\t$y2\n"; 
    my $xv = [$match_up_x1,$match_up_x2,$match_down_x2,$match_down_x1];
    my $yv = [$y1, $y1, $y2, $y2];
    my $points =$svg->get_path(
                     x=>$xv,y=>$yv,
                     -type=>'polyline',
                     -closed=>'true'
              );
    my $tag=$svg->polyline(
                     %$points,
                     style=>{
                        fill=>$color,
                        'fill-opacity' => 0.7
                     }
              );
}
}#end of function

sub plot_matches
{
my ($svg, $x1, $y1, $x2, $y2, $data) = @_;
my $scale = 500/5000;
#cut in sequence
my $upcut   = 0;
my $upcut_len= 0;
my $downcut = 0;
my $downcut_len = 0;
if (exists $data->{'up_cut'}){
   my @temp = split(',', $data->{'up_cut'});
   $upcut = $temp[0]; 
   $upcut_len = $temp[1];
}
if (exists $data->{'down_cut'}){
   my @temp = split(',', $data->{'down_cut'});
   $downcut = $temp[0];
   $downcut_len = $temp[1];
}


my @matches = split(';', $data->{'match'});
for (my $i=0; $i<@matches; $i++){
    print "$matches[$i]\n";
    my @array = split(',', $matches[$i]);
    my $match_up_x1 = $array[0];
    my $match_up_x2 = $array[1];
    my $match_down_x1 = $array[2];
    my $match_down_x2 = $array[3];
    my $match_strand = $array[4];
    print "Original matches: $match_up_x1, $match_up_x2, $match_down_x1, $match_down_x2\n";
    #convert position to relative cut point if need
    if ($upcut > 0){
       if ($match_up_x1 > $upcut){
           $match_up_x1 = ($match_up_x1 - $data->{'up_start'} - $upcut_len + 1)*$scale + $x1;
       }else{
           $match_up_x1 = ($match_up_x1 - $data->{'up_start'} + 1)*$scale + $x1;
       }
       if ($match_up_x2 > $upcut){
           $match_up_x2 = ($match_up_x2 - $data->{'up_start'} - $upcut_len + 1)*$scale + $x1;
       }else{
           $match_up_x2 = ($match_up_x2 - $data->{'up_start'} + 1)*$scale + $x1;       
       }
    }else{
       $match_up_x1 = ($match_up_x1 - $data->{'up_start'} + 1)*$scale + $x1;
       $match_up_x2 = ($match_up_x2 - $data->{'up_start'} + 1)*$scale + $x1; 
    }
 
    if ($downcut > 0){ 
       if ($match_down_x1 > $downcut){
           $match_down_x1 = ($match_down_x1 - $data->{'down_start'} - $downcut_len + 1)*$scale + $x2;
       }else{
           $match_down_x1 = ($match_down_x1 - $data->{'down_start'} + 1)*$scale + $x2;
       }
       if ($match_down_x2 > $downcut){
           $match_down_x2 = ($match_down_x2 - $data->{'down_start'} - $downcut_len + 1)*$scale + $x2;
       }else{
           $match_down_x2 = ($match_down_x2 - $data->{'down_start'} + 1)*$scale + $x2;        
       }
    }else{
       $match_down_x1 = ($match_down_x1 - $data->{'down_start'} + 1)*$scale + $x2;
       $match_down_x2 = ($match_down_x2 - $data->{'down_start'} + 1)*$scale + $x2;
    }

    # plot match
    print "match x: $match_up_x1\t$match_up_x2\t$match_down_x1\t$match_down_x2\n";
    print "match y: $y1\t$y1\t$y2\t$y2\n"; 
    my $xv = [$match_up_x1,$match_up_x2,$match_down_x2,$match_down_x1];
    my $yv = [$y1, $y1, $y2, $y2];
    my $points =$svg->get_path(
                     x=>$xv,y=>$yv,
                     -type=>'polyline',
                     -closed=>'true'
              );
    my $tag=$svg->polyline(
                     %$points,
                     style=>{
                        fill=>'gray'
                     }
              );
}

}

################################### sub for plot zoom in mping #######################
sub plot_zoom_in_mping
{
my ($svg, $x, $y, $pos, $data) =@_;
my $scale  = 500/5000;
my $start = 0;
my $end   = 0;
my $upcut   = 0;
my $upcut_len= 0;
my $downcut = 0;
my $downcut_len = 0;
my $zoom_in_y = $y;
my @zoom_in;
if ($pos =~ /up/){
   $start = $data->{'up_start'};
   $end   = $data->{'up_end'};
   if (exists $data->{'up_cut'}){
       my @temp = split(',', $data->{'up_cut'});
       $upcut = $temp[0]; 
       $upcut_len = $temp[1];
   }
   $zoom_in_y = $zoom_in_y - 45;
   @zoom_in = split(";", $data->{'up_zoom_in_mping'});
}elsif($pos =~ /down/){
   $start = $data->{'down_start'};
   $end   = $data->{'down_end'};
   if (exists $data->{'down_cut'}){
       my @temp = split(',', $data->{'down_cut'});
       $downcut = $temp[0];
       $downcut_len = $temp[1];
   }
   $zoom_in_y = $zoom_in_y - 45;
   @zoom_in = split(";", $data->{'down_zoom_in_mping'});
}

for (my $i=0; $i<@zoom_in; $i++){
   my @zoom_in_inf = split(',', $zoom_in[$i]);
   my $zoom_in_x = ($zoom_in_inf[0] - $start + 1)*$scale + $x;
   if ($pos =~ /up/ and $zoom_in_inf[0] > $upcut and $upcut > 0){
       $zoom_in_x = ($zoom_in_inf[0] - $start - $upcut_len + 1)*$scale + $x;
   }elsif($pos =~ /down/ and $zoom_in_inf[0] > $downcut and $downcut > 0){
       $zoom_in_x = ($zoom_in_inf[0] - $start - $downcut_len + 1)*$scale + $x;
   }
   my $zoom_in_left_seq  = $zoom_in_inf[2];
   my $zoom_in_right_seq = $zoom_in_inf[3];
   my $zoom_in_strand    = $zoom_in_inf[1];
   print "$zoom_in_inf[0]\t$start\t$x\t$zoom_in_x\t$zoom_in_y\n";
   plot_text($svg, $zoom_in_x-50, $zoom_in_y, 'middle', '10', $zoom_in_left_seq);
   plot_text($svg, $zoom_in_x+50, $zoom_in_y, 'middle', '10', $zoom_in_right_seq);
   plot_text($svg, $zoom_in_x, $zoom_in_y, 'middle', '10', 'mPing');
   plot_zoom_in_lines($svg, $zoom_in_x, $zoom_in_y)
}

}

################################### sub for plot_zoom in breakpoint signiature #######
sub plot_zoom_in
{
my ($svg, $x, $y, $pos, $data) =@_;
my $scale  = 500/5000;
my $start = 0;
my $end   = 0;
my $upcut   = 0;
my $upcut_len= 0;
my $downcut = 0;
my $downcut_len = 0;
my $zoom_in_y = $y;
my @zoom_in;
if ($pos =~ /up/){
   $start = $data->{'up_start'};
   $end   = $data->{'up_end'};
   if (exists $data->{'up_cut'}){
       my @temp = split(',', $data->{'up_cut'});
       $upcut = $temp[0];
       $upcut_len = $temp[1];
   }

   $zoom_in_y = $zoom_in_y - 45;
   @zoom_in = split(";", $data->{'up_zoom_in'});
}elsif($pos =~ /down/){
   $start = $data->{'down_start'};
   $end   = $data->{'down_end'};
   if (exists $data->{'down_cut'}){ 
       my @temp = split(',', $data->{'down_cut'});
       $downcut = $temp[0];
       $downcut_len = $temp[1];
   }
   $zoom_in_y = $zoom_in_y - 45;
   @zoom_in = split(";", $data->{'down_zoom_in'});
}
for (my $i=0; $i<@zoom_in; $i++){
   #print $zoom_in[$i],"\n";
   my @zoom_in_inf = split(',', $zoom_in[$i]);
   #print $zoom_in_inf[0],"\n";
   my $zoom_in_x = ($zoom_in_inf[0] - $start + 1)*$scale + $x;
   if ($pos =~ /up/ and $zoom_in_inf[0] > $upcut and $upcut > 0){
       $zoom_in_x = ($zoom_in_inf[0] - $start - $upcut_len + 1)*$scale + $x;
   }elsif($pos =~ /down/ and $zoom_in_inf[0] > $downcut and $downcut > 0){
       $zoom_in_x = ($zoom_in_inf[0] - $start - $downcut_len + 1)*$scale + $x; 
   }

   my $zoom_in_left_seq = $zoom_in_inf[1];
   my $zoom_in_right_seq = $zoom_in_inf[3];
   my $zoom_in_filler_seq = $zoom_in_inf[2];
   my $zoom_in_string = $zoom_in_left_seq.' '.$zoom_in_filler_seq.' '.$zoom_in_right_seq;
   print "$zoom_in_inf[0]\t$start\t$x\t$zoom_in_x\t$zoom_in_y\t$zoom_in_string\n"; 
   plot_text($svg, $zoom_in_x, $zoom_in_y, 'middle', '10', $zoom_in_string);
   plot_zoom_in_lines($svg, $zoom_in_x, $zoom_in_y)
} 
}

################################### sub for plot zoom in lines
sub plot_zoom_in_lines
{
my ($svg, $x, $y) = @_;
my $top_left = $x-80;
my $top_right= $x+80;
my $middle_y= $y + 30;
my $bottom_y= $y + 45; 
my $baseline=$svg->line(
     x1=>$x,y1=>$middle_y,
     x2=>$x,y2=>$bottom_y,
     style=>{stroke=>'black'}
   );
my $lefline=$svg->line(
     x1=>$x,y1=>$middle_y,
     x2=>$top_left,y2=>$y,
     style=>{stroke=>'black'}
   );   
my $rightline=$svg->line(
     x1=>$x,y1=>$middle_y,
     x2=>$top_right,y2=>$y,
     style=>{stroke=>'black'}
   );   

}

################################### convert position to relative cut point if need ####
#cut=cutposition,cutlength or 0,0
#up_start is the start position on chromosome
#match_up_x1 is the position on chromosome to convert
#scale=500/5000
#x1 is the start position in plot area
sub convert_position
{
my ($cut, $up_start, $match_up_x1, $scale, $x1) = @_;
my ($upcut, $upcut_len) = split(',', $cut);
if (1){
    if ($upcut > 0){
       if ($match_up_x1 > $upcut){
           $match_up_x1 = ($match_up_x1 - $up_start - $upcut_len + 1)*$scale + $x1;
       }else{
           $match_up_x1 = ($match_up_x1 - $up_start + 1)*$scale + $x1;
       }
    }else{
       $match_up_x1 = ($match_up_x1 - $up_start + 1)*$scale + $x1;
    }
}
return $match_up_x1
}


################################### sub for read ploting infromation###################
sub read_inf
{
my ($file) = @_;
my %hash;
my $rank=0;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/ or $_ =~/^#/);
    if ($_=~/^>(\d+)/){
        $rank = $1;
        next;
    }else{
        my @unit=split("\=", $_);
        #print $unit[0], "\t", $unit[1], "\n";
        $hash{$rank}{$unit[0]}=$unit[1];
    }
}
close IN;
return \%hash;
}


################################### sub for write svg to file##########################
sub writesvg {
my ($file, $svg)=@_;
open OUT, ">temp.svg" or die "can not open my file";
    print OUT $svg->xmlify;
close OUT;

my $left_arrow=<<DEF;
      <defs>
        <marker id="left_arrow" markerWidth="10" markerHeight="10" refX="0" refY="3" orient="180" markerUnits="strokeWidth" viewBox="0 0 10 10">
        <path d="M0,0 L0,6 L9,3 z" fill="black" />
        </marker>
      </defs>
DEF

my $right_arrow=<<DEF;
      <defs>
        <marker id="right_arrow" markerWidth="10" markerHeight="10" refX="0" refY="3" orient="0" markerUnits="strokeWidth" viewBox="0 0 10 10">
        <path d="M0,0 L0,6 L9,3 z" fill="black" />
        </marker>
      </defs>
DEF


open OUT, ">$file" or die "$!";
open IN, "temp.svg" or die "$!";
my $line = 0;
while(<IN>){
    chomp $_;
    next if ($_=~/^$/ or $_ =~/^#/);
    $line++;
    if ($line==3){
        print OUT "$_\n";
        print OUT "$left_arrow\n$right_arrow\n";
    }else{
        print OUT "$_\n";
    }
}
close IN;
system("/rhome/cjinfeng/BigData/software/draw/svg2xxx_release/svg2xxx $file -t pdf");
}
