// Global namespace
if (typeof site != 'undefined') {
    site._site = site;
}

//  Functions


// Simple JavaScript Templating
// John Resig - http://ejohn.org/ - MIT Licensed
(function(){
  var cache = {};

  this.tmpl = function tmpl(str, data){
    // Figure out if we're getting a template, or if we need to
    // load the template - and be sure to cache the result.
    var fn = !/\W/.test(str) ?
      cache[str] = cache[str] ||
        tmpl(document.getElementById(str).innerHTML) :

      // Generate a reusable function that will serve as a template
      // generator (and which will be cached).
      new Function("obj",
        "var p=[],print=function(){p.push.apply(p,arguments);};" +

        // Introduce the data as local variables using with(){}
        "with(obj){p.push('" +

        // Convert the template into pure JavaScript
        str
          .replace(/[\r\t\n]/g, " ")
          .split("<%").join("\t")
          .replace(/((^|%>)[^\t]*)'/g, "$1\r")
          .replace(/\t=(.*?)%>/g, "',$1,'")
          .split("\t").join("');")
          .split("%>").join("p.push('")
          .split("\r").join("\\'")
      + "');}return p.join('');");

    // Provide some basic currying to the user
    return data ? fn( data ) : fn;
  };
})();


// form input awaken
site.setInputConscious = function(t, off, on) {

    jQuery('input.' + t).each(function() {
        el = jQuery(this);
        tV = jQuery(this).val();

        jQuery(this).bind("click focus", function() {
            if (jQuery(this).val() === tV) {
                jQuery(this).val('');
                jQuery(this).toggleClass(on).toggleClass(off);
            }
        });

        jQuery(this).bind("focusout", function() {

            if (jQuery(this).val().length == 0) {
                jQuery(this).val(tV);
                jQuery(this).toggleClass(on).toggleClass(off);
            }

        });

        jQuery(this).closest('form').submit(function() {
            if (el.val() === tV) el.val('');
        });
    })
}

site.submitnewsletter = function(email) {
    $(email.messagecontainer).empty();  // in case of multiple clicks empty the message container.
    var subscribeToNewsletterUrl = "/subscribeToNewsletter";
    if (email.url != null){
        subscribeToNewsletterUrl = email.url+"/subscribeToNewsletter";
    }
    var subscriptionData = {"emailAddress" : email.address,
        "newsletters" : email.newsletters,
        "journalName" : null,
        "subscribed" : email.action};
    postnewsletters(subscribeToNewsletterUrl, subscriptionData, email);
}

site.submitnewsletters = function(email) {
    $(email.messagecontainer).empty();  // in case of multiple clicks empty the message container.
    var subscribeToNewsletterUrl = "/subscribeToNewsletter";
    if (email.url != null){
        subscribeToNewsletterUrl = email.url+"/subscribeToNewsletter";
    }
    var subscriptionData = {"emailAddress" : email.address,
        "newsletters" : email.newsletters,
        "journalName" : email.journalname,
        "subscribed" : email.action};
    postnewsletters(subscribeToNewsletterUrl, subscriptionData, email);
}

function postnewsletters(subscribeToNewsletterUrl, subscriptionData, email) {
    var data = $.ajax({
        url:subscribeToNewsletterUrl,
        type:"POST",
        data:subscriptionData,
        dataType:"json",
        timeout:5000,
        success:function (data, status) {
            if ((data.successful == true) && (email.action == true)) {
                $(email.messagecontainer).prepend('<p class="signup-success-msg success">Thank you - you have been added to the lists</p>');
            }
        },
        error:function (data, status, error) {
            if (data.status != 404) {
                var response = jQuery.parseJSON(data.responseText);

                if (typeof response.errors == "object" && response.errors != null && response.errors.toString().match("invalidEmail") != null && response.errors.toString().match("invalidEmail").length > 0) {
                    $(email.messagecontainer).prepend('<p class="error">Please enter a valid email address</p>');
                }
                else {
                    $(email.messagecontainer).prepend('<p class="error">There was a system error: ' + response.errors[0] + '</p>');
                }
            } else {
                $(email.messagecontainer).prepend('<p class="error">There was a system error: ' + data.status + " " + data.statusText + '</p>');
            }
        }
    });
}

site.videoPlaylist = function(args) {
    $f("flowplayer1", "http://releases.flowplayer.org/swf/flowplayer-3.2.7.swf", args).playlist("ul#clips1", {loop:true});
}

site.singleVideo = function(className) {
    flowplayer("a." + className, "http://releases.flowplayer.org/swf/flowplayer-3.2.7.swf", {

        // this is the player configuration.
        plugins:  {
            controls:  {
                volume: false
            }
        }
    });
}

site.toggleCheckBox = function() {
    jQuery('li.long').find('input:checked').each(function() {
        jQuery(this).parent('li').addClass('enabled');
    });
    jQuery('li.long').click(function() {
        jQuery(this).toggleClass('enabled');
        jQuery(this).find('input').each(function() {
            if (jQuery(this).is(':checked')) {
                jQuery(this).attr('checked', false);
            } else {
                jQuery(this).attr('checked', true);
            }
        })
    })
}

site.advancedSearchToggleCheck = function() {
    // full text long titles popup
    jQuery('.long label').each(function() {
        // capture abbreviated title in memory
        var abbreviatedText = jQuery(this).text();
        // replace abbreviated title with full title
        jQuery(this).mouseover(function() {
            jQuery(this).removeClass("clipped").addClass("full");
            jQuery(this).text($(this).attr('title'));
        });
        // replace full title with abbreviated title
        jQuery(this).mouseout(function() {
            jQuery(this).removeClass("full").addClass("clipped");
            jQuery(this).text(abbreviatedText);
        });
    });

        jQuery('.long a').click(function() {
        window.location.href = $(this).attr('href');
})

}

site.storeLink = function() {
    if (geoplugin_continentCode() == "EU"){
        $('#storeLink').html('<a href="http://biomedcentral.spreadshirt.net/?search_phrase=J'+ site.id + '">order merchandise on demand</a>');
    } else {
        $('#storeLink').html('<a href="http://biomedcentral.spreadshirt.com/?search_phrase=J'+ site.id + '">order merchandise on demand</a>');
        }
    };

// Intelligent pop-up resizing
site.resizeWindow = function() {

    // position popup in top left
    window.moveTo(0, 0);

    // init bodywidth + height
    var bodyWidth = 0;
    var bodyHeight = (document.body.clientHeight <= 80) ? jQuery(document).height() : document.body.clientHeight;
    //alert(bodyHeight);
    if (jQuery('#poparticle').width() > document.body.clientWidth) {
        bodyWidth = jQuery('#poparticle').width() + 50;
    }
    else {
        bodyWidth = document.body.clientWidth;
    }


// calculate toolbar height + width
    if (window.outerHeight) {
        var toolbarWidth = window.outerWidth - jQuery(window).width();
        var toolbarHeight = window.outerHeight - jQuery(window).height();

    }
    else {
        if (document.documentElement.clientWidth) {
            var clientW = document.documentElement.clientWidth;
            var clientH = document.documentElement.clientHeight;
            // resize the window, the new window will be smaller
            window.resizeTo(clientW, clientH);

            // calculate the difference between the client and the total size in the new window
            var frameW = clientW - document.documentElement.clientWidth;
            var frameH = clientH - document.documentElement.clientHeight;

            var toolbarWidth = frameW
            var toolbarHeight = frameH;
        }
    }
//alert(toolbarHeight);


// calculate total height and width
    var totalHeight = bodyHeight + toolbarHeight;
    var totalWidth = bodyWidth + toolbarWidth;
    //alert(totalHeight);

    var finalWidth = 0;
    var finalHeight = 0;

    // if totalwidth is bigger than screen, minimise to screenwidth
    if (totalWidth > screen.availWidth || bodyWidth > screen.availWidth) {
        finalWidth = screen.availWidth;
        ////console.log('body width=' + bodyWidth + ' final width=' + finalWidth + ' horiz-max\n');
    }
    // otherwise wrap content to screen width
    else {
        finalWidth = totalWidth;
        ////console.log('body width=' + bodyWidth + ' final width=' + finalWidth + ' horiz-min\n');
    }

    // if body + toolbars is bigger than screen, window should maximise to screen.
    if (totalHeight > screen.availHeight || bodyHeight > screen.availHeight) {
        // maximise
        finalHeight = screen.availHeight;
        ////console.log('body height=' + bodyHeight + ' final height=' + finalHeight + ' vert-max');

    }
    // else innerwindow should be just enough to show content
    else {
        finalHeight = totalHeight;
        ////console.log('body height=' + bodyHeight + ' final height=' + finalHeight + ' vert-min');
    }

    window.resizeTo(finalWidth, finalHeight);

}

// will collapse an accordion if anchor on the h2 matches hash in url
site.collapseAccordionByAnchor = function() {
    if (window.location.hash.length > 0) {
        if(site.page == "article_stats") {
            if(window.location.hash == "#citations-biomedcentral") {
                $('h2#citations').click();
                $('h3#citations-biomedcentral').click();        
            } else {
            $('h2' + window.location.hash).click();
        }
    }
    }
}

function format(selected, total) {
    return selected + '/' + total + ' selected';
}

site.initializeAccordion = function(options) {
    defaultOptions = {
        activeClass: 'new-action',
        speed:'fast'/*'slow', 'fast', 'def'*/,
        defaultActive:0,
        showDefault:true
    }

    this.options = jQuery.extend(defaultOptions, options);

// jQuery('.wrap-in', '.wrap').find('table').wrap('<div/>');
    jQuery('h2', '.wrap').removeClass('active');
    jQuery('h2', '.wrap:has(table.' + this.options.activeClass + ')').addClass('active');
    jQuery('.wrap-in', '.wrap').find('table.' + this.options.activeClass).parent('div').prev('h3').addClass('active');
    jQuery('.wrap-in', '.wrap').find('table').not('.' + this.options.activeClass).parent('div').prev('h3').removeClass('active'); // removed .hide() after parent(div)
    jQuery('.wrap-in', '.wrap').not(':has(table.' + this.options.activeClass + ')').hide();
    jQuery('.wrap-in').each(
            function(index, element) {
                var parent = jQuery(element).parent();
                var note = parent.find('.note').show();
            });

    var self = this;

    jQuery('.wrap').find('h2').click(
            function() {
                jQuery(this).parent().find('div').
                        eq(0).slideToggle(self.options.speed,
                        function() {
                            jQuery(this).prev('h2').toggleClass('active')
                        })
            }
            )
    jQuery('.wrap-in', '.wrap').find('h3').click(
            function() {
                jQuery(this).next('div').slideToggle(self.options.speed,
                        function() {
                            jQuery(this).prev('h3').toggleClass('active')
                        })
            }
            )

    if (this.options.showDefault) {
        if (jQuery('.wrap:has(table.' + this.options.activeClass + ')').length == 0) {
            jQuery('.wrap').eq(this.options.defaultActive).find('h2').addClass('active');
            jQuery('.wrap').eq(this.options.defaultActive).find('.wrap-in').show();
        }
    }

    if(this.options.hideAll === true) {
        jQuery('h2').removeClass('active');
        jQuery('.wrap-in').hide();
    }
}


site.accordionHover = function() {
    $('.accordion .wrap-in h3').hover(function() {
        $(this).toggleClass('hover');
    });
}


site.ieSelectExpand = function(el, original, expand) {
    var option = el + ' option';
    if (jQuery.browser.msie) {

        jQuery(el).focusin(function(event) {
            $(this).width(expand);
            $(this).blur();
        });
        jQuery(el).change(function() {
            $(this).width(original);
        });
    }
}

/* generates tooltip on all alt tags within el - must also include tooltip.js on styles.vm for this
 to work */
function tooltip(el) {
    jQuery(el).tooltip({
        track: true,
        delay: 0,
        showURL: false,
        fixPNG: true,
        showBody: " - ",
        top: -20,
        left: 15
    });
}

/* used to toggle browse by subject accordion */
site.toggler = function(el) {
    var toggler = jQuery(el);

    toggler.each(function() {
    if($(this).parent('li').attr('id') != window.location.hash.replace('#','')) {
        $(this).add(jQuery('i', $(this)).removeClass('active'));
        $(this).next().hide();
    } else {
        $(this).addClass('active');
        $(this).next().show();
    }
    });



    toggler.click(function() {
        jQuery(this).add(jQuery('i', this)).toggleClass('active');
        jQuery(this).next().slideToggle();
        jQuery(this).toggleClass('selected');
    }
            );

    if (site.section == "journals" || (site.section == "about" && site.subSection == "editors")) {
        toggler.mouseover(function(e) {
            jQuery(this).css('background-color', '#F1F1F1');
        });


        toggler.mouseout(function(e) {
            if ((jQuery(this).hasClass('selected')) != true) {
                jQuery(this).css('background-color', 'transparent');
            }
        });

        tooltip('.tooltip');



    }
}

site.addTextAreaMaxLengthListener = function() {

    $('textarea[maxlength]').bind('paste cut keyup input blur load',
            function(e) {
                //get the limit from maxlength attribute
                var limit = parseInt($(this).attr('maxlength'));
                //get the current text inside the textarea
                var text = $(this).val();
                //count the number of characters in the text
                var chars = text.length;

                var label = $("label[for='" + this.id + "']");
                //var label = $('label[for="userInterests.interests"]');

                if (label) {
                    var new_limit = limit - chars;
                    label.text(new_limit);
                }
                // do not merge with if above
                if (chars >= limit && $('span.error').text().length === 0) {
                    $('.textarea-count').show();
                    // disable form submission whilst char > limit
                    jQuery('#update-interests').attr('disabled', 'true').attr('style', 'background-color: #ccc');
                    jQuery('input').bind('keydown keypress keyup', function(event) {
                        if (event.which == 13) {
                            event.preventDefault();
                        }
                    });
                }
                else {
                    $('.textarea-count').hide();
                    jQuery('#update-interests').removeAttr('disabled').removeAttr('style');
                }
            }
            );
}

// collapsing list utility
// toggle visibility of element that is a sibling of an element with class .collapser
site.toggleList = function() {
    if (jQuery('.collapser').siblings('ul').attr('class') != 'hidebeforeload') {
        jQuery('.collapser').siblings('ul').toggle(); // show list if js is disabled
    }
    jQuery('.collapser').click(function() {
        jQuery(this).siblings('ul').toggle();
        jQuery(this).children('i').toggleClass('active');

    });
}

// Remove advert block when there is no ad inside. Show 'advert' text when there is an ad
site.hideAdverts = function() {
    jQuery('.google-ad,.skyscraper,.advert').each(function() {

        // if there is no ad hide container
        if (jQuery(this).find('img').length > 0 && jQuery(this).find('img').attr('src').match('empty.gif')) {
            if (jQuery.browser.msie && jQuery.browser.version <= 7) {
                jQuery(this).css({'visibility':'hidden'});
            }
            else if (jQuery.browser.msie && jQuery.browser.version == 8) {
                jQuery(this).css({'visibility':'hidden', 'height':'0'});
            }

            else {
                jQuery(this).css({'display':'none', 'height':'0'});
            }
        } // otherwise show container
        else if (jQuery(this).children('dd').height() > 20) {
            jQuery(this).children('dt').css('display', 'block');
        }
        // display advertisement bg image for top ad
        if (jQuery(this).hasClass('wide') && jQuery(this).children('dd').height() > 20) { // for wide ads which us bg image
            jQuery(this).addClass('advert-text');
        }
    });
}

site.abtesting = function() {

    if(typeof ABalytics == "object") {
        ABalytics.applyHtml();
        var pos = $(".ab_test1").text() || "bottom_center";
        $('#bmc_cancer_mpu dl.google-ad dd a').detach().prependTo(pos + " dl.google-ad dd"); 
        $(pos).show(); 
    }
}

// If there is an element with id "scroll-here" than it will be focused
site.scrollToTop = function() {
    if (jQuery('#scroll-here').length == 1) {
        window.location = "#scroll-here";
    }
}

function popup(url, n, w, h) {
    var popup_window = 0;

    if (!popup_window.top) {
        popup_window = window.open(url, n, 'width=' + w + ',height=' + h + ',resizable=yes,scrollbars=yes,toolbar=no,location=no,directories=no,status=no,menubar=no,copyhistory=no');
    }
    popup_window.focus();
}

site.addParameterIfSkipCache = function() {
    if (/skipCache/.test(window.location.href)) {
        return "?skipCache=true";
    } else {
        return "";
    }
}

site.getParameterIfSkipCache = function() {
    if (/skipCache/.test(window.location.href)) {
        return "&skipCache=true";
    } else {
        return "";
    }
}

site.reverse = function(str) {
    if (!str) return '';
    var revstr = '';
    for (i = str.length - 1; i >= 0; i--)
        revstr += str.charAt(i)
    return revstr;
}
//  onLoadJson(id of panel in which to insert feed, number of results to show, url to fetch, panel )

site.onLoadJson = function (feedUri, length, panel, heading, template, prefix, args) {

    // template to use, falls back to default if not specified
    var jsTemplate = template || "feed_template";
    var jsPrefix = prefix || "/webapi/1.0/blogs/";
    var webApiUrl;

    // using regex to match word "search", to treat searchuri's differently
    var containsSearch = feedUri.match(/\bsearch/g);
    var twitterMatch = feedUri.match(/\btwitter/g);

    var mostviewed = feedUri.match(/\bmost_viewed/g);

    if (containsSearch !== null && twitterMatch === null) {
        webApiUrl = site.abbreviation + feedUri + site.addParameterIfSkipCache();
    }
    else if (mostviewed !== null) {
        webApiUrl = site.abbreviation + "/webapi/1.0/" + feedUri + ".json" + site.addParameterIfSkipCache();
    }
    else {
        var jsPrefix = prefix || "/webapi/1.0/blogs/";
        webApiUrl = site.abbreviation + prefix + length.toString() + '?feedUrl=' + encodeURIComponent(site.reverse(feedUri)) + site.getParameterIfSkipCache();
    }


    // parsing json object
    function displayFeed(panel, data) {
        data.heading = heading;

        // data.params is a container object for any optional parameters passed to onLoadJson
        data.params = args || {}; // add extra paramaters if passed to onLoadJson

        // extend params with quoteBySearchParams if they are available in this call
        if (typeof quoteBySearchParams == "object") {
            jQuery.extend(data.params, quoteBySearchParams);
        }

        jQuery('.' + panel + '_feed').html(tmpl(jsTemplate, data));
        // console.log(data);
    }

    // config parameters
    var options = {
        type: "GET",
        url: webApiUrl,
        dataType: "json",
        timeout: 60000, // 6 seconds
        success: function(data, textStatus) {
            if ((typeof data.entries == "object" && data.entries.length > 0) || (typeof data.articles == "object")) {    // only display if there is something to show
                displayFeed(panel, data);
            } else if (typeof console == "object") {
                console.log(data.title + " is hidden because it returned no results. It's an " + typeof data.entries);
            } // otherwise show what there is for debugging
        },
        error: function(xhr, textStatus, errorThrown) {
            //console.log(webApiUrl+ ' returned error: ' +errorThrown + '\n and status:' + textStatus + ' ' + xhr.responseText);
        },
        complete: function(xhr, textStatus) {
            //console.log(webApiUrl+ '\n complete - returned status: '+' ' +textStatus);
        }
    };

    // fetch json
    jQuery.ajax(options);


}

//json volume handler
site.jsonVolumesEventHandler = function() {
    jQuery('.jsonVolumeLoader').unbind();
    jQuery('.jsonVolumeLoader').bind('click', function(event) {

        var volumeObj = jQuery(this);
        var volume = jQuery(this).attr('id');
        event.preventDefault();

        var webApiUrl = site.abbreviation + "/browseVolume.json?volume=" + volume + site.getParameterIfSkipCache();

        var options = {
            type: "GET",
            url: webApiUrl,
            dataType: "json",
            timeout: 60000, // 6 seconds
            beforeSend: function (XMLHttpRequest) {
                // startAjaxSpinner();
            },
            success : function(data, textStatus) {
                var volumeNumber = data.articleCollectionDataList[0].volumeNumber;
                updateVolumeList(volumeNumber, data);
                jQuery("#" + volumeNumber).children('h3').addClass("active");
                jQuery('div.volume_' + volumeNumber).show();
                //stopAjaxSpinner();
            },
            error : function (XMLHttpRequest, textStatus, errorThrown) {
                // jQuery('.ajax-error').slideDown("fast").fadeOut(6000);
            },
            complete: function (XMLHttpRequest, textStatus) {
                //stopAjaxSpinner();
            }
        };
        //if (jQuery.trim(jQuery("div.volume_" +  volume).html()) != "") {
        if (jQuery(this).children('h3').hasClass("active")) {
            jQuery("div.volume_" + volume).html("");
            jQuery(this).children('h3').removeClass("active");

        }
        else {
            jQuery.ajax(options);

        }
    });


    function startAjaxSpinner() {
        jQuery('.ajax-spinner').html("<img src=\"/images/icons/ajax-loader.gif\" alt=\"\" />").show();
    }

    function stopAjaxSpinner() {
        jQuery('.ajax-spinner').hide();
    }

    function updateVolumeList(vol, data) {
        jQuery("div.volume_" + vol).html(tmpl("render_issue_volume", data));
    }

}

 site.recordOutboundLink = function(link, category, action) {
    _gat._getTrackerByName()._trackEvent(category, action);
    setTimeout('document.location = "' + link.href + '"', 100);
  }

//json archive handler
site.jsonArchiveEventHandler = function() {
    jQuery('.jsonArchiveLoader').unbind();
    jQuery('.jsonArchiveLoader').bind('click', function(event) {
        var volumeNumber = jQuery(this).attr('id');
        event.preventDefault();

        var webApiUrl = site.abbreviation + "/browseVolume.json?volume=" + volumeNumber + site.getParameterIfSkipCache();
        var options = {
            type: "GET",
            url: webApiUrl,
            dataType: "json",
            timeout: 15000, // 6 seconds
            beforeSend: function (XMLHttpRequest) {
                // startAjaxSpinner();
            },
            success : function(data, textStatus) {
                // var yearOfPublicationNumber = data.year;
                updateVolumeList(volumeNumber, data);
                jQuery("#" + volumeNumber).children('h3').addClass("active");
                jQuery('div.archive_' + volumeNumber).show();
                //stopAjaxSpinner();
            },
            error : function (XMLHttpRequest, textStatus, errorThrown) {
                // jQuery('.ajax-error').slideDown("fast").fadeOut(6000);
            },
            complete: function (XMLHttpRequest, textStatus) {
                //stopAjaxSpinner();
            }
        };

        if (jQuery.trim(jQuery("div.archive_" + volumeNumber).html()) != "") {
            jQuery("div.archive_" + volumeNumber).html("");
            jQuery(this).children('h3').removeClass("active");
        }
        else {
            jQuery.ajax(options);
        }
    });


    function startAjaxSpinner() {
        jQuery('.ajax-spinner').html("<img src=\"/images/icons/ajax-loader.gif\" alt=\"\" />").show();
    }

    function stopAjaxSpinner() {
        jQuery('.ajax-spinner').hide();
    }

    function updateVolumeList(volumeNumber, data) {
        jQuery("div.archive_" + volumeNumber).html(tmpl("render_archive_volume", data));
        if(site.page == "archive") { jQuery("div.archive_" + volumeNumber).find('.tooltip').tooltip();  }
    }

}

//json volume decade handler
site.jsonDecadeEventHandler = function() {
    jQuery('.jsonIssueDecadeLoader').unbind();
    jQuery('.jsonIssueDecadeLoader').bind('click', function(event) {
        var decade = jQuery(this).attr('id');
        event.preventDefault();

        var webApiUrl = site.abbreviation + "/browseDecade.json?decadeStart=" + decade + site.getParameterIfSkipCache();

        var options = {
            type: "GET",
            url: webApiUrl,
            dataType: "json",
            timeout: 15000, // 6 seconds
            beforeSend: function (XMLHttpRequest) {
                // startAjaxSpinner();
            },
            success : function(data, textStatus) {
                var decadeNumber = Math.floor(data.volumes[0].startYear / 10) + "0";
                updateDecadeList(decadeNumber, data);
                site.jsonVolumesEventHandler();
                jsonArchiveEventHandler();
                jQuery("#" + decadeNumber).children('h2').addClass("active");
                //stopAjaxSpinner();
            },
            error : function (XMLHttpRequest, textStatus, errorThrown) {
                // jQuery('.ajax-error').slideDown("fast").fadeOut(6000);
            },
            complete: function (XMLHttpRequest, textStatus) {
                //stopAjaxSpinner();
            }
        };

        if (jQuery.trim(jQuery("div.decade_" + decade).html()) != "") {
            jQuery("div.decade_" + decade).html("");
            jQuery(this).children('h2').removeClass("active");
        }
        else {
            jQuery.ajax(options);
        }
    });


    function startAjaxSpinner() {
        jQuery('.ajax-spinner').html("<img src=\"/images/icons/ajax-loader.gif\" alt=\"\" />").show();
    }

    function stopAjaxSpinner() {
        jQuery('.ajax-spinner').hide();
    }

    function updateDecadeList(decade, data) {
        jQuery("div.decade_" + decade).html(tmpl("render_decade", data));
    }

}

site.jsonEventHandler = function(preload) {
    var cache = {};
    jQuery('.jsonloader').bind('click', function(event) {
        var preCacheKey = jQuery('#article-pane-nav').find('.current').children().attr('id');
        if (!cache[preCacheKey]) {
            cache[preCacheKey] = jQuery('#articles').html();
        }
        linkId = jQuery(this).attr('id');
        var parameters = jQuery(this).attr('jsonParams');
        event.preventDefault();
        var webApiUrl = site.abbreviation + "/webapi/1.0/" + linkId + ".json" + site.addParameterIfSkipCache();
        if (parameters) {
            webApiUrl += "?" + parameters;
        }
        var options = {
            type: "GET",
            url: webApiUrl,
            dataType: "json",
            timeout: 6000, // 6 seconds
            success : function(data, textStatus) {
                update(linkId, data, null);
            },
            error : function (XMLHttpRequest, textStatus, errorThrown) {
                if(preload == false) {
                    jQuery('.ajax-error').slideDown("fast").fadeOut(6000);
                }
            }
        };
        if (cache[linkId]) {
            update(linkId, null, cache[linkId]);
        } else {
            jQuery.ajax(options);
        }
    });

    function update(container, data, html) {
         if(preload == false) {    
            displayRss();
            updateArticleList(container, data, html);
            makeTabCurrent(container);
        }
    }

    function displayRss() {
        jQuery('ul#article-pane-nav li').find('a.rss').hide();
        jQuery('#' + linkId).siblings().children('a.rss').show();
    }

    function makeTabCurrent(linkId) {
        jQuery('.article-panel-tabs').removeClass('current');
        jQuery('#' + linkId).parent().addClass('current');
    }

    function updateArticleList(linkId, data, html) {
        if (!html) {
            cache[linkId] = tmpl(linkId + "_template", data);
        }
        jQuery('#articles').html(cache[linkId]);
        return html;
    }

    if(preload == true) {
		if(jQuery('.jsonloader').length > 1) {
			/* this was not working in safari and android browsers
			for(var i = 1; i < jQuery('.jsonloader').length; i++) {
			   jQuery('.jsonloader')[i].click(); 
			}
			*/
			$('#article-pane-nav a.jsonloader').each( function(){
				$(this).click()
			});
       }
    }

}


// Hide all rss icons apart from first
site.hideRssIcons = function() {
    jQuery('ul#article-pane-nav li').eq(0).siblings().find('a.rss').hide();
}

// Advanced Search
jQuery(function() {
    var optionsToggler = jQuery('p.options a.options-toggler');

    optionsToggler.click(
            function() {
                jQuery(this).add(jQuery("i", this)).toggleClass('active');
                jQuery(this).toggleClass('active');
            });
});

// IE 6 and 7 input:focus
jQuery(document).ready(function() {
    if ( $.browser.msie ) {

        jQuery("input[type='text']:not(.placeholder)").focus(function() {
            jQuery("input:[text]").css('color', '#000000');
        });
        jQuery("div.tertiary-content fieldset.block-form .text").focus(function() {
            jQuery(this).css('color', '#000');
        });

        jQuery("div.tertiary-content fieldset.block-form .text").blur(function() {
            jQuery(this).css('color', '#D1D1D1');
        });

    }

});



site.initialiseScrollable = function() {


    // if .random class present gallery items will be shuffled
    //--------------------------------------------------------
    var gall = jQuery('div#thumbs div');
    if (jQuery("div.scrollable").hasClass("random")) {

        // grab and randomise order of gallery items
        
        gall.sort(function() {
            return 0.5 - Math.random()
        });
        gall.sort(function() {
            return 0.5 - Math.random()
        });

        // remove existing unshuffled items
        jQuery('div#thumbs div').remove();

        // add shuffled items
        for (var i = 0; i < gall.length; i++) {
            jQuery('div#thumbs').prepend(gall[i]);
        }

        // display gallery
        jQuery("div.scrollable").css("visibility", "visible");

    }


    // Setting the interval (default is 10 seconds)
    var scrollableInterval = 10000;
    if($("div.scrollable").hasClass("slow")) {
        scrollableInterval = 5000;
    }

    // initialize scrollable
    jQuery("div.scrollable").scrollable({
        size: 1,
        items: '#thumbs',
        hoverClass: 'hover',
        loop: (gall.length > 1) ? true : false,
        circular: (gall.length > 1) ? true : false
    }).navigator({
                    // show dots
                    navi: "div.navi"
                 }).autoscroll({
                                autoplay: (!site.plugins.mobileTopSection) ? true : false,
                                interval: scrollableInterval
                               });

    // single image widget should not show arrows
    // navigation is usually hardcoded so this cannot be an option, must manually detach the elements.
    if(gall.length == 1) {
        $("div.scrollable-container div.navi, div.scrollable-container a.browse").detach();
    }

    // arrow show/hide behaviour
    $('div.scrollable-container').bind('mouseover', function() {
		//if(!site.plugins.mobileTopSection){
			$(this).children('a.left, a.right, div.navi').removeClass('disabled');
		//}
    }).bind('mouseout', function() {
		//if(!site.plugins.mobileTopSection) {
			$(this).children('a.left, a.right, div.navi').addClass('disabled');
		//}
		/*
		else {
			$(this).children('a.left, a.right').addClass('disabled');
		}
		*/
	});
	
	//running for mobile to centrally align navi buttons by setting width
	if(site.plugins.mobileTopSection) {
		var naviButtons = $("div.navi a");
		var scrollableCount = naviButtons.size();
			//console.log(scrollableCount);
		var getNaviWidth = naviButtons.outerWidth(true);
			//console.log(getNaviWidth);
		var newNaviWidth = scrollableCount * getNaviWidth+"px";
			//console.log(newNaviWidth);
		$("div.navi").css("width", newNaviWidth);
	}
}


site.setUpPlayer = function(){
    $("a.flowplayer").each(function(){
        $f(this.id, $(this).attr("flash"));
    });
}

site.updateTwitterPanel = function() {
    getTwitters('twitter_updates', {
        id: 'biomedcentral',
        count: 3,
        ignoreReplies: true,
        newwindow: true,
        clearContents: false,
        timeout: 10,
        onTimeoutCancel: true,
        callback: function() {
            jQuery("#twitter_panel").css('display', 'block');
        }
    });
};

// Andrew Urquhart : CountDownPro Timer : www.andrewu.co.uk/clj/countdown/pro/
function CD_M(strTagId) {
    var objMeta = document.getElementsByTagName("meta");
    if (objMeta && objMeta.length) {
        for (var i = 0; i < objMeta.length; ++i) {
            if (objMeta.item(i).scheme == strTagId) {
                var name = objMeta.item(i).name;
                var content = objMeta.item(i).content;
                if (name.indexOf("mindigits") > 0 || name.indexOf("hidezero") > 0) {
                    window[strTagId][name] = parseInt(content, 10)
                } else {
                    window[strTagId][name] = content
                }
            }
        }
    }
}
;
function CD_UD(strContent, objW) {
    objW.node.innerHTML = strContent
}
;
function CD_T(strTagId) {
    var objNow = new Date();
    var objW = window[strTagId];
    if (objW.msoffset) {
        objNow.setMilliseconds(objNow.getMilliseconds() + objW.msoffset)
    }
    ;
    CD_C(objNow, objW);
    if (objW.intEvntDte <= objNow.valueOf()) {
        if (objW.event_functionhandler && typeof window[objW.event_functionhandler] == "function") {
            window[objW.event_functionhandler](new Date(objW.intEvntDte));
            objW.event_functionhandler = ""
        }
        if (objW.event_msg || objW.event_redirecturl) {
            var msg = "<span id=\"" + strTagId + "_complete\">" + objW.event_msg + "</span>";
            if (objW.event_redirecturl) {
                location.href = objW.event_redirecturl
            } else if (objW.event_audio_src) {
                var strMimeType = objW.event_audio_mimetype;
                var audioObject = "<object style=\"visibility:hidden;\" id=\"MediaPlayer\" width=\"2\" height=\"2\" data=\"" + objW.event_audio_src + "\" type=\"" + strMimeType + "\"></object>";
                CD_UD(msg + audioObject, objW)
            } else {
                CD_UD(msg, objW)
            }
            ;
            return
        }
    }
    ;
    setTimeout("if(typeof CD_T=='function'){CD_T(\"" + strTagId + "\")}", 1100 - objNow.getMilliseconds())
}
;
function CD_C(objNow, objW) {
    var intMS = objW.intEvntDte - objNow.valueOf();
    if (intMS <= 0) {
        intMS *= -1
    }
    ;
    var intD = Math.floor(intMS / 864E5);
    intMS = intMS - (intD * 864E5);
    var intH = Math.floor(intMS / 36E5);
    intMS = intMS - (intH * 36E5);
    var intM = Math.floor(intMS / 6E4);
    intMS = intMS - (intM * 6E4);
    var intS = Math.floor(intMS / 1E3);
    var strTmp = CD_F(intD, "d", objW) + CD_F(intH, "h", objW) + CD_F(intM, "m", objW) + CD_F(intS, "s", objW);
    CD_UD(strTmp, objW)
}
;
function CD_F(intData, strPrefix, objW) {
    if (intData == 0 && objW[strPrefix + "_hidezero"]) {
        return ""
    }
    ;
    var strResult = "" + intData;
    var intMinDigits = objW[strPrefix + "_mindigits"];
    if (intData.toString().length < intMinDigits) {
        strResult = "0000000000" + strResult;
        strResult = strResult.substring(strResult.length, strResult.length - intMinDigits)
    }
    if (intData != 1) {
        strResult += objW[strPrefix + "_units"]
    } else {
        strResult += objW[strPrefix + "_unit"]
    }
    ;
    return objW[strPrefix + "_before"] + strResult + objW[strPrefix + "_after"]
}
;
function CD_Parse(strDate) {
    var objReDte = /(\d{4})\-(\d{1,2})\-(\d{1,2})\s+(\d{1,2}):(\d{1,2}):(\d{0,2})\s+GMT([+\-])(\d{1,2}):?(\d{1,2})?/;
    if (strDate.match(objReDte)) {
        var d = new Date(0);
        d.setUTCFullYear(+RegExp.$1, +RegExp.$2 - 1, +RegExp.$3);
        d.setUTCHours(+RegExp.$4, +RegExp.$5, +RegExp.$6);
        var tzs = (RegExp.$7 == "-" ? -1 : 1);
        var tzh = +RegExp.$8;
        var tzm = +RegExp.$9;
        if (tzh) {
            d.setUTCHours(d.getUTCHours() - tzh * tzs)
        }
        if (tzm) {
            d.setUTCMinutes(d.getUTCMinutes() - tzm * tzs)
        }
        ;
        return d
    } else {
        return NaN
    }
}
;
function CD_Init() {
    var strTagPrefix = "countdown";
    var objElem = true;
    if (document.getElementById) {
        for (var i = 1; objElem; ++i) {
            var strTagId = strTagPrefix + i;
            objElem = document.getElementById(strTagId);
            if (objElem && (typeof objElem.innerHTML) != 'undefined') {
                var strDate = objElem.innerHTML;
                var objDate = CD_Parse(strDate);
                if (!isNaN(objDate)) {
                    var objW = window[strTagId] = new Object();
                    objW.intEvntDte = objDate.valueOf();
                    objW.node = objElem;
                    objW.servertime = "";
                    objW.d_mindigits = 1;
                    objW.d_unit = " day";
                    objW.d_units = " days";
                    objW.d_before = "";
                    objW.d_after = " ";
                    objW.d_hidezero = 0;
                    objW.h_mindigits = 2;
                    objW.h_unit = "h";
                    objW.h_units = "h";
                    objW.h_before = "";
                    objW.h_after = " ";
                    objW.h_hidezero = 0;
                    objW.m_mindigits = 2;
                    objW.m_unit = "m";
                    objW.m_units = "m";
                    objW.m_before = "";
                    objW.m_after = " ";
                    objW.m_hidezero = 0;
                    objW.s_mindigits = 2;
                    objW.s_unit = "s";
                    objW.s_units = "s";
                    objW.s_before = "";
                    objW.s_after = " ";
                    objW.s_hidezero = 0;
                    objW.event_msg = "";
                    objW.event_audio_src = "";
                    objW.event_audio_mimetype = "";
                    objW.event_redirecturl = "";
                    objW.event_functionhandler = "";
                    CD_M(strTagId);
                    if (objW.servertime) {
                        var objSrvrTm = CD_Parse(objW.servertime);
                        if (isNaN(objSrvrTm)) {
                            objElem.innerHTML = strDate + "**";
                            continue
                        } else {
                            objW.msoffset = parseInt((objSrvrTm.valueOf() - (new Date()).valueOf()) / 1000, 10) * 1000
                        }
                    } else {
                        objW.msoffset = 0
                    }
                    ;
                    CD_T(strTagId);
                    if (objElem.style) {
                        objElem.style.visibility = "visible"
                    }
                } else {
                    objElem.innerHTML = strDate + "<a href=\"http://andrewu.co.uk/clj/countdown/pro/\" title=\"CountdownPro Error:Invalid date format used,check documentation (see link)\">*</a>"
                }
            }
        }
    }
}
if (window.attachEvent) {
    window.attachEvent('onload', CD_Init)
} else if (window.addEventListener) {
    window.addEventListener("load", CD_Init, false)
} else {
    window.onload = CD_Init
}
;


// detect zoom event - http://mlntn.com/2008/12/11/javascript-jquery-zoom-event-plugin/
jQuery.fn.zoom = function(fn) {
    // Set handler for keyboard zooming in Firefox, IE, Opera, Safari.
    // This is the only valid case of browser-specific code I've ever seen -JM
    jQuery(document).keydown(function(e) {
        switch (true) {
            case jQuery.browser.mozilla || jQuery.browser.msie :
                if (e.ctrlKey && (
                        e.which == 187 || // =/+ (zoom in [FF, IE])
                                e.which == 189 || // -   (zoom out [FF, IE])
                                e.which == 107 || // +   (numpad) (zoom in [FF, IE])
                                e.which == 109 || // -   (numpad) (zoom out [FF, IE])
                                e.which == 96 || // 0   (reset zoom [FF, IE])
                                e.which == 48     // 0   (numpad) (reset zoom [IE, FF, Opera])
                        )) fn();
                break;
            case jQuery.browser.opera :
                // Opera requires CTRL to be pressed for reset (using num 0)
                if (
                        e.which == 43 || // +   (numpad) (zoom in [Opera, Safari])
                                e.which == 45 || // -   (numpad) (zoom out [Opera, Safari])
                                e.which == 42 || // *   (numpad) (reset zoom [Opera])
                                (e.ctrlKey && e.which == 48) // 0   (numpad) (reset zoom [FF, IE, Opera])
                        ) fn();
                break;
            case jQuery.browser.safari :
                // Use e.metaKey for the Apple key
                if (e.metaKey && (
                        e.charCode == 43 || // +   (numpad) (zoom in [Opera, Safari])
                                e.charCode == 45    // -   (numpad) (zoom out [Opera, Safari])
                        )) fn();
                break;
        }
        return;
    });

    // Set handler for scrollwheel zooming in IE
    jQuery(document).bind('mousewheel', function(e) {
        if (e.ctrlKey) fn();
    });

    // Set handler for scrollwheel zooming in Firefox
    jQuery(document).bind('DOMMouseScroll', function(e) {
        if (e.ctrlKey) fn();
    });
};

site.quickSearchBoxResize = function() {

    // ascertain free width
    var input = jQuery('input#searchTerms').width();
    var text = jQuery('fieldset.search span').width();
    var newwidth = 355 - text;

    // add extra width to input panel
    jQuery('.branding-inner input#searchTerms').width(newwidth + 'px');
}

site.uploadFileField = function() {
    // used to simulate file upload input field

    $('input.upload-field').wrap('<div class="upload-container right">');
    $('.upload-container').before('<input class="file text"/>');

    $('input#file').change(function() {
        $('input.file.text').attr('value', $('input#file').attr('value').replace('C:\\fakepath\\', ''))
    });

    $('input.file.text').focus(function() {
        $('input#file').click();
    });

    $('input.upload-field').css('display', 'block');
}

/* function for handling drop down flash ads which need to overlay content */
site.popOverBannerAd = function(popOverAd) {
    if ($.browser.msie && ($.browser.version == "7.0" || $.browser.version == "6.0")) {
        $(popOverAd).mouseover(function() {
            $(this).css({'overflow' : 'visible'});
        })

        $('body').mousemove(function(e) {
            if (e.pageY > 300 || e.pageX < 250 || e.pageX > 1000) {
                $(popOverAd).css({'overflow' : 'hidden'});
            }
        });
    }
    else if ($.browser.mozilla) {
        $('#branding').prepend('<div id="ff-hitarea"></div>');

        $('#ff-hitarea').mouseover(
                function() {
                    $(popOverAd).addClass('popover');
                    $(popOverAd).parent('dl').css('height', '90px');
                }).mouseout(function() {
            $(popOverAd).removeClass('popover');
        });
    }
    else {
        $(popOverAd).mouseover(
                function() {
                    $(this).addClass('popover');
                    $(this).parent('dl').css('height', '90px');
                }).mouseout(function() {
            $(this).removeClass('popover');
        });

    }

}

site.toggleList = function() {
    if (jQuery('.collapser').siblings('ul').attr('class') != 'hidebeforeload') {
        jQuery('.collapser').siblings('ul').toggle(); // show list if js is disabled
    }
    jQuery('.collapser').click(function() {
        jQuery(this).siblings('ul').toggle();
        jQuery(this).children('i').toggleClass('active');

        if($("span.toggletext").length > 0) {  jQuery(this).children("a.collapser span").toggle();   }

        if (site.page == "microsite") {
            jQuery('#showmore').css('display', 'none');

            /* gateways signup box  */
            if($("#signup-developingworld").length > 0) {
                $("#signup-developingworld").submit(function(e) {
                    e.preventDefault();
                    site.submitnewsletter({address: $(this).find("input#email").attr("value"), newsletters: null, action:true, messagecontainer: "#message-box"} );
                });
            }
        }
    });
}

site.updatesPanel = function() {
    jQuery('#twitter_updates_tab').bind('click', function() {
        jQuery(this).addClass('current');
        jQuery('#fromtheblog_tab').removeClass('current');
        jQuery('#twitter_updates').css('display', 'block');
        jQuery('#fromtheblog').css('display', 'none');
    });

    jQuery('#fromtheblog_tab').bind('click', function() {
        jQuery(this).addClass('current');
        jQuery('#twitter_updates_tab').removeClass('current');
        jQuery('#twitter_updates').css('display', 'none');
        jQuery('#fromtheblog').css('display', 'block');
    });
}

    site.citationArticlesAjax = function() {
        //ajax for citation articles
        var ajaxArticleDiv = jQuery("#citations_biomedcentral_articles");
        ajaxArticleDiv.hide();

        jQuery("#about-citedon-biomed-link").click(function() {
            jQuery(this).parent().click();
            return false;
        });

        jQuery("#citations-biomedcentral").css("cursor", "pointer");
        jQuery("#citations-biomedcentral").click(function() {
            if (ajaxArticleDiv.hasClass("article_loaded")) {
                ajaxArticleDiv.toggle();
                ajaxArticleDiv.css('display','block');
                return false;
            }
            var url = $("#urlForAjax").html();
            jQuery.ajax({
                url: url,
                timeout:300000,
                beforeSend:function() {
                    ajaxArticleDiv.html('<img src="/images/loading.gif" alt="loading"/>');
                    ajaxArticleDiv.show();
                },
                error:function(XMLHttpRequest, textStatus, errorThrown) {
                    ajaxArticleDiv.html('Sorry, we cannot load the articles at the moment. Please try later.');
                    //console.log(errorThrown);
                },
                success: function(data) {
                    ajaxArticleDiv.html(data);
                    ajaxArticleDiv.addClass("article_loaded");
                    if(jQuery('a.summary').length > 0) {
                           jQuery('a.summary').click(function(){
                            var summaryToggled=jQuery(this).toggleClass('active').parents('.nav').eq(0).next('div.summary-toggled').eq(0);
                            var self=this;
                            summaryToggled.slideToggle('fast');
                            return false;
                        });
                }
                }
            });
        });
        
        
        
    }
// on load

function loadInstitutionRecognitionDetails() {    
    jQuery.ajax({
        url:"/newjournal/institutionRecognitionDetails",
        dataType:"html",
        success:function (data) {
            if (data != null && data.match("institutional_recognition_panel") != null) {
                $("#welcome_text_without_institution_panel").hide();
                $("#welcome-message-panel").html(data);
                if(site.id==10283) { stripPortalLinks(); }
            }
        },
        error:function (request, errorType, errorThrown) {
            // $("#welcome-message-panel").html('<h4 class="error">Sorry, an error occurred: ' + errorType + '</h4>')
        },
        cache:false
    });
}

function loadInstitutionRecognitionHTMLFragment() {

    // if cookie isnt present force load it by calling a non-varnished resource
    if(site.allCookies.getItem("bmccookie") == null) {
        $.ajax({url:"/logon"})
    }

	var name = "institution_id";
	var cookie = site.allCookies.getItem("bmccookie");
	var getinstid = getValueFromURLParam(name, cookie);
	var path = '/inst/panel/'+getinstid;
	if(getinstid !== null && getinstid !=="" ) {
		getInstPanel(path);
	} else {
		//showing generic template for now.
		loadInstitutionRecognitionDetails();
	}
}

function getInstPanel(requestUrl) {
	jQuery.ajax({
		url: requestUrl,
		dataType:"html",
		success:function (data) {
            if (data !== null && data.match("institutional_recognition_panel") !== null) {
				$("#welcome_text_without_institution_panel").hide();
				$("#welcome-message-panel").html(data);
                if(site.id==10283) { stripPortalLinks(); }
			}
		},
		error:function (request, errorType, errorThrown) {
            $("#welcome-message-panel").html('<h4 class="error">Sorry, an error occurred: ' + errorType + '</h4>');
		},
		cache:false
	});
}

function getValueFromURLParam(name, string) {
	name = name.replace(/[\[]/,"\\\[").replace(/[\]]/,"\\\]");
	var regexS = "[\\?&]"+name+"=([^&#]*)";
	var regex = new RegExp( regexS );
	var results = regex.exec( string );
	if(results == null)
		return "";
	else
	return results[1];
}

function stripPortalLinks() {
    // needs to be abstracted for all journal welcome panels (currently only springerplus)
    if($("#welcome-message-panel").length > 0) {
        $("#welcome-message-panel ul.intro li:contains('SpringerOpen journals')").children("a").attr("href","http://springeropen.com/journals");
        var inst = $("#welcome-message-panel a[href*='inst']").attr("href");
        $("#welcome-message-panel a[href*='inst']").attr("href","http://springeropen.com/"+inst);
    }
}

jQuery(document).ready(function() {
    if (site.portal != true) { // don't resize SO searchbox
        site.quickSearchBoxResize();
    }

    if (site.page == "microsite") {
        site.updatesPanel();
    }

    site.scrollToTop();
    site.toggleList();
    site.hideRssIcons();
    site.addTextAreaMaxLengthListener();
    // ajax
    if (site.section == "home") { // run on homepage only
		if((site.portal || site.id == 10283) && window.bmcIsMobile != "mobile") { 
			// load the institution recognition details

			//DEV-3296
			
            if(site.id == 10283) {
			 loadInstitutionRecognitionHTMLFragment();
            }
			else {
			loadInstitutionRecognitionDetails();
		}
		}
	
        if (quoteSearchUri) {

            if (typeof heading == "undefined") {
                site.onLoadJson(quoteSearchUri, 0, 'quote_by_search', '', 'quotes_by_search_template', '/webapi/1.0/search/');
            }
            else if (site.id == 3001 || site.id == 3008 || site.id == 10184) { // petsko's column and musings on GM
                site.onLoadJson(quoteSearchUri, 0, 'quote_by_search', heading, 'custom_quote_by_search_column', '/webapi/1.0/search/');
            }
            else {
                site.onLoadJson(quoteSearchUri, 0, 'quote_by_search', heading, 'quotes_by_search_template', '/webapi/1.0/search/');
            }
        }

        site.ieSelectExpand('select#journalDropList', '217px', '370px');

    }


    if (site.plugins.scrollable) {
        site.initialiseScrollable();
    }
    
    if (site.section == "journals" || site.section == "about") {
        site.toggler('.toggler');
    }

    if (site.page == "indexing" || site.page == "article_stats" || site.page == "funderpoliciespage" || site.page == "institutionalpolicies") {
        if(site.page == "article_stats") {
            site.initializeAccordion({showDefault: true});
            site.citationArticlesAjax ();
        } else {
        site.initializeAccordion();
    }
    }

    if (site.page == "ordinst_form" || site.page == "free_trials" || site.page == "advancedsearch" || site.page == "quick_reg_contact") {
        site.toggleCheckBox();
    }

    if (site.page == "advancedsearch") {
        site.advancedSearchToggleCheck();
    }

    if (site.page == "biography") {
        site.uploadFileField();
    }
    
    if (site.page == "advocacy") {
        site.storeLink();
    }   

    if(site.page == "home" && site.id == 2019) {
        site.abtesting();
    }

    // first time the ajax requests are loaded in background and cached
    site.jsonEventHandler(true);
    // second time time the click event is attached to tabs
    site.jsonEventHandler(false);

    site.jsonVolumesEventHandler();
    site.jsonArchiveEventHandler();
    site.jsonDecadeEventHandler();
    // function needs to be modified for archive page
    if (site.page != "archive") {
        site.accordionHover();
    } 
    
    if(site.page == "archive") { $('.tooltip').tooltip(); }

	/* for journals only
		$("#signupForUpdates").submit(function(e) {
			e.preventDefault();
			site.submitnewsletter({address: $(this).find("input#email").attr("value"), newsletters: "JOURNAL_SPECIFIC_NEWSLETTER", action:true, messagecontainer: "#message-box"} );
		}); 
		}
	*/

    /*  sign up for cases-updates about/updates  */
    if($("#signupForUpdates").length > 0 && site.id == "10121") {
        $("#signupForUpdates").submit(function(e) {
            e.preventDefault();
            site.submitnewsletter({address: $(this).find("input#email").attr("value"), newsletters: "cases-updates", action:true, messagecontainer: "#message-box"} );
        });
    }
	
    /**
       check if dbstats search_cookei parameter is available from default_dbstats_layout.vm
    */
    if(typeof search_cookie!= 'undefined')
    {
        //  readingDBStatsCookie() sets the value of "request_Id" in  default_dbstats_layout.vm
        readingDBStatsCookie();

            if(typeof request_Id != 'undefined')
            {
                    var par1="bytes sent via SQL*Net to client";
                    var par2="physical read bytes";
                    var par3="consistent gets";
                     jQuery.ajax(
                             {
                                     url : "/newjournal/dbstats?content_type=json&stats_name=" + par1 +"&stats_name="+ par2 +"&stats_name="+ par3+ "&request_id="+request_Id,
                                     success : function(data)
                                     {
                                        var url = $('#more_dbstats_link a').attr('href');
                                         for(j=0; j < data.length; j++)
                                         {
                                             $.each(data[j], function(key,value)
                                                     {
                                                             if(value.indexOf(par3)!= -1){
                                                                 $("#consistent_gets").text(getDBStats(value));
                                                             }
                                                             if(value.indexOf(par2)!= -1){
                                                                 $("#physical_reads").text(getDBStats(value));
                                                             }
                                                             if(value.indexOf(par1)!= -1){
                                                                 $("#SQL_NET").text(getDBStats(value));
                                                             }
                                                     }
                                             );
                                         }
                                    }
                             }
                     );

                     function getDBStats(o)
                     {
                         var value1 = JSON.stringify(o);
                         var splitString = value1.split(",");
                         var i = splitString.length;
                         var temp = splitString[i-1];
                         return(temp.replace(/[^a-zA-Z 0-9]+/g,''))
                     }
            }
    }
});
// webkit cannot handle document.ready with width() or height()
jQuery(window).load(function() {
    site.hideAdverts();

    // if there is a collapsing banner ad
    if (parseInt($('.google-ad.wide #OAS_RMF_Top_FLASH').attr('height')) > 90) {
        site.popOverBannerAd(('.google-ad.wide dd'));
    }

    // resize pop-up windows (chrome has issue with resizing and doc.load)
    if (jQuery("body").attr('id') == "popup") {
        site.resizeWindow();
    }
})

site.loadPanelByJson = function(configs) {
    for (var i = 0; i < configs.length; i++) {
        site.onLoadJson(configs[i].feedUri, configs[i].feedLength, configs[i].panel, configs[i].heading, configs[i].template, configs[i].prefix, configs[i].args);
    }
}

site.renderCollapsable = function(header) {
    var isActive = header.hasClass('active');
    var collapseBody = header.next(".collapse-body");
    if (isActive) {
        collapseBody.show();
    } else {
        collapseBody.hide();
    }
}
site.initializeUserProfilePage = function() {
    jQuery('.collapse-header').click(function() {
        var collapsableHeader = $(this);
        collapsableHeader.toggleClass('active');
        site.renderCollapsable(collapsableHeader);
    });

    $('.collapse-header').each(function(index) {
        site.renderCollapsable($(this));
    });
}

function toggleContactDetails(contactStatus, contactDetails) {
    var checked = $("input[name=" + contactStatus + "]").is(':checked');
    $("input[name=" + contactDetails + "]").attr('disabled', !checked);
}

/*
function submitManuscript(ac) {
    jou_id = document.getElementById('journalId').value;
    if (ac == 'submit' && jou_id == '10280') document.submit_manuscript.action = 'mailto:comc@connected.ly';
    else if (ac == 'instructions' && jou_id == '10280') document.submit_manuscript.action = "#";
}*/

function submitManuscriptGateway() {
    var siteurl = "http://www." + $('#submit_manuscript_gateway option:selected').val() + "/manuscript"
    if ($('#submit_manuscript_gateway option:selected').val() !="") {
        window.location.href=siteurl;
    }
}

function submitInstructionsGateway() {
    var siteurl = "http://www." + $('#submit_manuscript_gateway option:selected').val() + "/authors/instructions"
    if ($('#submit_manuscript_gateway option:selected').val() !="") {
        window.location.href=siteurl;
    }
}

function submitform() {
    document.submit_manuscript.action = '/authors/instructions/#submitManuscript';
    var currentElement = document.createElement("input");
    currentElement.setAttribute("type", "hidden");
    currentElement.setAttribute("name", "actionPerformed");
    currentElement.setAttribute("value", "viewScope");
    document.getElementById('submit_manuscript').appendChild(currentElement);
    document.getElementById('submit_manuscript').submit();
}

function LoadInParent(link) {

    $('li' + link).css('background-color','#f2f2f2');
       window.location = link;
    window.focus();
}

function setIframeHeight(iframe) {
    if (iframe) {
        var iframeWin = iframe.contentWindow;
        if (iframeWin.document.body) {
            iframe.height = iframeWin.document.documentElement.scrollHeight+10;
        }
    }
}

function endsWith(str, suffix) {
    return str.indexOf(suffix, str.length - suffix.length) !== -1;
}

site.initializeJumBrowseForm = function() {
    $('#jumpBrowse').submit(function() {
        if ($('#jumpToVolume').attr('value') == '') {
            $('#showMandatoryVolumeErrorMessage').show();
            return false;
        }
    });
}


site.initializeInterestsForm = function() {
    $('h2.active').removeClass('active').siblings('.wrap-in').css('display', 'none');
    $('h2.expanded').addClass('active').siblings('.wrap-in').css('display', 'block');
    $('h2.active').removeClass('expanded');
}

site.initializeLightbox = function() {
    /* Fancy box for IP /libraries/free_trials*/
 $("a#whatismyip").fancybox();

    /*
     *   Group Images and Navigation within Lightbox
     */

    $("a[rel=photos_group]").fancybox({
        'transitionIn'        : 'none',
        'transitionOut'        : 'none',
        'titlePosition'     : 'over',
        'titleFormat'        : function(title, currentArray, currentIndex, currentOpts) {
            return '<span id="fancybox-title-over">Image ' + (currentIndex + 1) + ' / ' + currentArray.length + (title.length ? ' &nbsp; ' + title : '') + '</span>';
        }
    });

    ////console.log('lightbox has run');
}



site.loadJournalSection = function() {
    $("#journalSectionId").change(function() {
        jQuery.ajax({
            url : "/author/manuscript/details/sectionDescription",
            dataType: "html",
            data : {
                "journalSectionId" : $("#journalSectionId").val()
            },
            success : function(data) {
                $("#sectiondescription_result").show();
                $("#sectiondescription_result").html(data);
                $("#subsectiondescription_result").hide();

                },
            error : function(request, errorType, errorThrown) {
                $("#sectiondescription_result").html('<h4 class="error">Sorry, an error occurred: ' + errorType + '</h4>');
            },
            cache : false,
            timeout : 60000
        });
    });

    $("#sectiondescription_result").show();
    $("#loadSubsections").hide();
    $("#scope").hide();
    $("#viewScope").hide();
}

site.loadJournalSubSectionDescription = function() {
    $("#journalSubSectionId").change(function() {
        jQuery.ajax({
            url : "/author/manuscript/details/subSectionDescription",
            dataType: "html",
            data : {
                "journalSubSectionId" : $("#journalSubSectionId").val()
            },
            success : function(data) {
                $("#subsectiondescription_result").html(data);
                $("#subsectiondescription_result").show();
                },
            error : function(request, errorType, errorThrown) {
                $("#subsectiondescription_result").html('<h4 class="error">Sorry, an error occurred: ' + errorType + '</h4>');
            },
            cache : false,
            timeout : 60000
        });
    });

    $("#subsectiondescription_result").show();
    $("#loadSubsections").hide();
    $("#scope").hide();
    $("#viewSubsectionScope").hide();
}

site.loadJournalSubSection = function() {
    $("#journalSectionId").change(function() {
        jQuery.ajax({
            url : "/author/manuscript/details/subSections",
            dataType: "html",
            data : {
                "journalSectionId" : $("#journalSectionId").val()
            },
            success : function(data) {
                $("#journalSubSectionId").html(data);
                if  ($('#journalSubSectionId option').size() <= 1) {
                    $("#journalSubSectionId").hide();
                    if (!site.loadSubjects_shouldDisableSubjects()) {
                        site.loadSubjects_enableSubjects();
                    }
                } else {
                    $("#journalSubSectionId").show();
                    site.loadSubjects_disableSubjects();
                }
            },
            error : function(request, errorType, errorThrown) {
                alert("Error while loading subsections")
            },
            cache : false,
            timeout : 60000
        });
    });
    if  ($('#journalSubSectionId option').size() <= 1) {
        $("#journalSubSectionId").hide();
    }

     $("#loadSubsections").hide();
     $("#scope").hide();
}

site.loadSubjects_doLoadSubjects = function(journalId, journalSectionId) {
    $.getJSON("/author/subjects/primary?journalId=" + journalId + "&journalSectionId=" + journalSectionId, function(data) {
        var items = [];
        $.each(data, function(key, val) {
            items.push('<option value="' + key + '">' + val + '</option>');
        });
        items.sort(function(a, b) {
            var o1 = $(a).text();
            var o2 = $(b).text();
            return o1 < o2 ? -1 : o1 > o2 ? 1 : 0;
        });
        $('select[name="primarySubject"] option[value!="0"]').remove();
        $('select[name="primarySubject"]').append(items.join(''));

        var primarySubject = $('input:hidden[name=primarySubject]').val();
        $('select[name=primarySubject]').val(primarySubject);
        // set selected option, if one is chosen
        if (site.elementExists('input[name=primarySubject]')) {
            var primarySubjectId = $('input[name=primarySubject]').val();
            var selector = "select[name=primarySubject] option[value=" + primarySubjectId + "]";
            $(selector).attr('selected', true);
        }
    });
    $.getJSON("/author/subjects/secondary?journalId=" + journalId + "&journalSectionId=" + journalSectionId, function(data) {
        var items = [];
        $.each(data, function(key, val) {
            items.push('<option value="' + key + '">' + val + '</option>');
        });
        items.sort(function(a, b) {
            var o1 = $(a).text();
            var o2 = $(b).text();
            return o1 < o2 ? -1 : o1 > o2 ? 1 : 0;
        });
        $('select[name^="secondarySubject"] option[value!="0"]').remove();
        $('select[name^="secondarySubject"]').append(items.join(''));


        var secondarySubjectOne = $('input:hidden[name=secondarySubjectOne]').val();
        $('select[name=secondarySubjectOne]').val(secondarySubjectOne);

        var secondarySubjectTwo = $('input:hidden[name=secondarySubjectTwo]').val();
        $('select[name=secondarySubjectTwo]').val(secondarySubjectTwo);

        // set selected option, if one is chosen
        if (site.elementExists('input[name=secondarySubjectOne]')) {
            var secondarySubjectIdOne = $('input[name=secondarySubjectOne]').val();
            var selector = "select[name=secondarySubjectOne] option[value=" + secondarySubjectIdOne + "]";
            $(selector).attr('selected', true);
        }

        if (site.elementExists('input[name=secondarySubjectTwo]')) {
            var secondarySubjectIdTwo = $('input[name=secondarySubjectTwo]').val();
            var selector = "select[name=secondarySubjectTwo] option[value=" + secondarySubjectIdTwo + "]";
            $(selector).attr('selected', true);
        }
    });


}

site.loadSubjects_shouldLoadSubjects = function() {
    if ($('select[name="journalSectionId"]').length) {
        if ($('select[name="journalSubSectionId"]').length && (parseInt($('select[name="journalSubSectionId"]').val()) > 0)) {
            return true;
        }
        if ($('select[name="journalSubSectionId"]').length == 0 && $('select[name="journalSectionId"]').length && (parseInt($('select[name="journalSectionId"]').val()) > 0)) {
            return true;
        }
    }
    return site.loadSubjects_getJournalSectionId() > 0;
}

site.loadSubjects_shouldDisableSubjects = function() {
    if (site.elementExists('select#journalSectionId') && !site.elementHasValue('select[name="journalSectionId"]')) {
        return true;
    }
    if (site.elementExists('select#journalSubSectionId:visible') && !site.elementHasValue('#journalSubSectionId:visible')) {
        return true;
    }
    return false;
}

site.elementExists = function(selector) {
    return $(selector).length > 0
}

site.elementHasValue = function(selector) {
    return parseInt($(selector).val()) > 0;
}

site.loadSubjects_disableSubjects = function() {
    $('select[name="primarySubject"]').attr('disabled', 'disabled');
    $('select[name="secondarySubjectOne"]').attr('disabled', 'disabled');
    $('select[name="secondarySubjectTwo"]').attr('disabled', 'disabled');
}

site.loadSubjects_enableSubjects = function() {
    $('select[name="primarySubject"]').removeAttr('disabled');
    $('select[name="secondarySubjectOne"]').removeAttr('disabled');
    $('select[name="secondarySubjectTwo"]').removeAttr('disabled');
}

site.loadSubjects_getJournalSectionId = function() {
    if ($('#journalSubSectionId').hasClass('submitted')) {
        return parseInt($('#journalSubSectionId').attr('title'));
    }
    if (site.elementExists('#journalSubSectionId:visible')) {
        return parseInt($('#journalSubSectionId').val());
    }
    return parseInt($('#journalSectionId').val());
}

site.loadSubjects_getJournalId = function() {
    return $('input[name="journalId"]').attr('value') || "";
}

site.loadSubjects = function(){
    if (site.loadSubjects_shouldDisableSubjects()) {
        site.loadSubjects_disableSubjects();
    }

    if (site.loadSubjects_shouldLoadSubjects()) {
        site.loadSubjects_doLoadSubjects(site.loadSubjects_getJournalId(), site.loadSubjects_getJournalSectionId());
    }
    $('#journalSectionId').change(function(event) {
        if (!site.loadSubjects_shouldDisableSubjects()) {
            site.loadSubjects_enableSubjects();
        } else {
            site.loadSubjects_disableSubjects();
        }
        site.loadSubjects_doLoadSubjects(site.loadSubjects_getJournalId(), site.loadSubjects_getJournalSectionId());
    });
    $('#journalSubSectionId').change(function(event) {
        if (!site.loadSubjects_shouldDisableSubjects()) {
            site.loadSubjects_enableSubjects();
        } else {
            site.loadSubjects_disableSubjects();
        }
        site.loadSubjects_doLoadSubjects(site.loadSubjects_getJournalId(), site.loadSubjects_getJournalSectionId());
    });
}

/* FAQ and browse/articles scroll up/down to refactor as a function (greg) */
jQuery('a.options-toggler').click(function() {
    jQuery('.checkbox-toggle').toggle();
    jQuery(this).toggleClass('active').parents('.options').eq(0).next().slideToggle({complete:function() {
        jQuery('.display-controls').toggle();
    }});
    return false;
});

/* browse by subject anchor toggler */
site.toggleAnchorAccordion = function() {
    if(jQuery(location.hash).length > 0) {
        $(location.hash).children('h2').toggleClass('active').siblings('ul:hidden').show();
    }
}

site.youtubelightbox = function(video) {

    jQuery(video.container).click(function(e) {
        e.preventDefault();
        jQuery.fancybox({
            'padding'       : 0,
            'autoScale'     : false,
            'transitionIn'  : 'none',
            'transitionOut' : 'none',
            'title'         : this.title,
            'width'     : video.width,
            'height'        : video.height,
            'href'          : this.href.replace(new RegExp("watch\\?v=", "i"), 'v/'),
            'type'          : 'swf',
            'swf'           : {
                'wmode'     : 'transparent',
                'allowfullscreen'   : 'true'
            }
        });

        return false;
    });
}

site.widen = function(panel,size) {
    $(panel).css("width",size);
}

/* cookie getter/setter */
site.allCookies = {
          getItem: function (sKey) {
            if (!sKey || !this.hasItem(sKey)) { return null; }
            return unescape(document.cookie.replace(new RegExp("(?:^|.*;\\s*)" + escape(sKey).replace(/[\-\.\+\*]/g, "\\$&") + "\\s*\\=\\s*((?:[^;](?!;))*[^;]?).*"), "$1"));
          },
          setItem: function (sKey, sValue, vEnd, sPath, sDomain, bSecure) {
            if (!sKey || /^(?:expires|max\-age|path|domain|secure)$/i.test(sKey)) { return; }
            var sExpires = "";
            if (vEnd) {
              switch (vEnd.constructor) {
                case Number:
                  sExpires = vEnd === Infinity ? "; expires=Tue, 19 Jan 2038 03:14:07 GMT" : "; max-age=" + vEnd;
                  break;
                case String:
                  sExpires = "; expires=" + vEnd;
                  break;
                case Date:
                  sExpires = "; expires=" + vEnd.toGMTString();
                  break;
              }
            }
            document.cookie = escape(sKey) + "=" + escape(sValue) + sExpires + (sDomain ? "; domain=" + sDomain : "") + (sPath ? "; path=" + sPath : "") + (bSecure ? "; secure" : "");
          },
          removeItem: function (sKey, sPath) {
            if (!sKey || !this.hasItem(sKey)) { return; }
            document.cookie = escape(sKey) + "=; expires=Thu, 01 Jan 1970 00:00:00 GMT" + (sPath ? "; path=" + sPath : "");
          },
          hasItem: function (sKey) {
            return (new RegExp("(?:^|;\\s*)" + escape(sKey).replace(/[\-\.\+\*]/g, "\\$&") + "\\s*\\=")).test(document.cookie);
          },
          keys: /* optional method: you can safely remove it! */ function () {
            var aKeys = document.cookie.replace(/((?:^|\s*;)[^\=]+)(?=;|$)|^\s*|\s*(?:\=[^;]*)?(?:\1|$)/g, "").split(/\s*(?:\=[^;]*)?;\s*/);
            for (var nIdx = 0; nIdx < aKeys.length; nIdx++) { aKeys[nIdx] = unescape(aKeys[nIdx]); }
            return aKeys;
          }
        };

// Prototype augmentations 
String.prototype.capitalCase = function() {
    return this.charAt(0).toUpperCase() + this.slice(1);
}

//check if cookie are enabled
    function are_cookies_enabled()
    {
    	var cookieEnabled = (navigator.cookieEnabled) ? true : false;

    	if (typeof navigator.cookieEnabled == "undefined" && !cookieEnabled)
    	{
    		document.cookie="testcookie";
    		cookieEnabled = (document.cookie.indexOf("testcookie") != -1) ? true : false;
    	}
    	return (cookieEnabled);
    }

jQuery(document).ready(function() {

    if(site.page == "collections") {
        if(jQuery('a.summary').length > 0) {
               jQuery('a.summary').click(function(){
                var summaryToggled=jQuery(this).toggleClass('active').parents('.nav').eq(0).next('div.summary-toggled').eq(0);
                var self=this;
                summaryToggled.slideToggle('fast');
                return false;
            });
        }
    }

     if(site.page == "collections" && jQuery("a#youtube").length > 0) {
        site.youtubelightbox({container: "a#youtube", width: 680, height:495})
     }

    if(jQuery("#journal-launch-notification").length > 0) {
        jQuery("#journal-launch-notification").submit(function(e) {
            e.preventDefault();
            site.submitnewsletter({address: $(this).find("input#email").attr("value"), newsletters: "newjournals", action:true, messagecontainer: "#message-box"} );
        });
    }

        if(jQuery("#journal-launch-notification-springer").length > 0) {
        jQuery("#journal-launch-notification-springer").submit(function(e) {
            e.preventDefault();
            site.submitnewsletter({address: $(this).find("input#email").attr("value"), newsletters: "new-so-journals", action:true, messagecontainer: "#message-box"} );
        });
    }


    if(jQuery("#articleSignupForUpdates").length > 0) {
        jQuery("#articleSignupForUpdates").submit(function(e) {
            e.preventDefault();
            site.submitnewsletters({address: $(this).find("input#email").attr("value"), newsletters: "bmcupdates, JOURNAL_SPECIFIC_NEWSLETTER", journalname: $(this).find("input#journalName").attr("value"), action:true, messagecontainer: "#message-box"} );
        });
    }
	
		// close sign up box on article pages
	if (site.section == "browse" ){
		$('.close').click(function(){
			$('.article-alert-signup-div').hide();
		})
	}
	
    // editor/manuscript pages are too narrow, right column shifts below so expand container
    if($("table.admin_content").length > 0) {
        site.widen("#content","1275px");
    }

    site.loadPanelByJson(jsonLoaders.list());
    
    jQuery("#journalListId").css('width', '334px');
    if (site.page == "emailPreferences" || site.page == "details") {
        site.initializeInterestsForm();    //for profile/my_interests_form.vm and user_interests_form.vm
    }
    site.initializeJumBrowseForm();
    toggleContactDetails('globalContactStatus', 'globalContactDetails');
    toggleContactDetails('journalContactStatus', 'journalContactDetails');
    site.initializeUserProfilePage();
    site.setInputConscious('placeholder', 'slept', 'awaken');
    if (site.page == "indexing" || site.page == "article_stats") {
        site.collapseAccordionByAnchor();
    }

    if(site.section == "my" && site.page == "myManuscripts") {
        //change subsection link (this will appear if the editor has suggested a section when reject and transfer)
        $('.changeSection').click(function(){
                    //hide section name
                    $(this).parent().children(".sectionName").hide();

                    //hide view scope link
                    $(this).parent().children('.showSectionDescription').hide();

                    //hide the change link
                    $(this).hide();

                    //show section drop dawn
                    var select = $(this).parent().children('.journalSectionId');
                    select.show();

                    //hide the subsection  tr
                    var subsectionTr = $(this).parent().parent().siblings(".subsection");
                    subsectionTr.hide();

                    //populate the hidden field - section
                    var op1 = $(select).children("option:selected"), id = $(op1).attr("class"), optionValue1 =  $(op1).attr("value");
                    $("input#articleSection-" + id).attr("value",optionValue1);

                    //populate the hidden field - subSection
                    $("input#articleSubsection-" + id).attr("value","");

                });

        // Load Journal Section description;
        $('.journalSectionId').change(function() {
            var select = this;
            var scopelink = $(select).siblings("span").parent().children('.showSectionDescription');
            var desc =  $(select).siblings("span").children('.sectiondescription_result') || "";
            var subsecdesc = $(select).siblings("span").children('.subsectiondescription_result') || "";
            $(scopelink).children('a').removeClass('viewScope') ;

            if(typeof String.prototype.trim !== 'function') {
              String.prototype.trim = function() {
                return this.replace(/^\s+|\s+$/g, '');
              }
            }
            if($(select).val() != "0") {

                jQuery.ajax({
                    url : "/author/manuscript/details/simpleSectionDescription",
                    dataType: "html",
                    data : {
                        "journalSectionId" : $(select).val()
                    },
                    success : function(data) {
                        $(desc).html(data);
                        $(desc).hide();
                        $(subsecdesc).hide();
                        if($(desc).text().trim().length > 0) {
                            $(scopelink).show();
                        } else {
                            $(scopelink).hide();
                        }
                        },
                    error : function(request, errorType, errorThrown) {
                        $(desc).html('<h4 class="error">Sorry, an error occurred: ' + errorType + '</h4>');
                    },
                    cache : false,
                    timeout : 60000
                });
            } else {
                $(scopelink).hide();
            }
        });


        // Load Journal SubSection
        $('.journalSectionId').change(function() {
            var select = this;
            var subsection =  $(this).parents("tr").siblings(".subsection");
            var subsectionselect = $(this).parents("tr").siblings(".subsection").find("select");
            var subsectionName = $(subsectionselect).parent().children(".subSectionName");
            var subsectionChangeLink = $(subsectionselect).parent().children(".changeSubsection");
            var subsedescccontainer = $(subsection).find(".showSubSectionDescription");

            //hide subsection name and change button if present
            subsectionName.hide();
            subsectionChangeLink.hide();

            //populate the hidden field - section
            var op1 = $(select).children("option:selected"), id = $(op1).attr("class"), optionValue1 =  $(op1).attr("value");
            $("input#articleSection-" + id).attr("value",optionValue1);

            //hide section error
            $("#articleSectionError-" + id).hide();

            //populate the hidden field - subSection
            var op2 = $(subsectionselect).children("option:selected"), optionValue2 =  $(op2).attr("value");

            $(this).siblings(".showSectionDescription a").removeClass("viewScope");
            $(subsedescccontainer).hide();

            if($(select).val() != "0") {
                jQuery.ajax({
                    url : "/author/manuscript/details/subSections",
                    dataType: "html",
                    data : {
                        "journalSectionId" : $(select).val()
                    },
                    success : function(data) {
                        $(subsectionselect).html(data);
                        if  ($(subsectionselect) && $(subsectionselect).children("option").length <= 1) {
                            $(subsection).hide();
                            $("input#articleSubsection-" + id).attr("value","");
                        } else {
                            $(subsection).show();
                            $(subsectionselect).show();
                            $("input#articleSubsection-" + id).attr("value","0");
                        }
                    },
                    error : function(request, errorType, errorThrown) {
                        alert("Error while loading subsections")
                    },
                    cache : false,
                    timeout : 60000
                });
            } else {
                $(subsection).hide();
            }
        });

        // Load Journal SubSectionDescription
        $('.journalSubSectionId').bind("change", function() {
            var select = this;
            var subSectionOption = $(this).children("option:selected");

            var scopelink = $(select).siblings("span").parent().children(".showSubSectionDescription");

            $(scopelink).children('a').removeClass('viewScope') ;
            var desc =  $(select).siblings("span").children('.subsectiondescription_result') || "";

            //populate the hidden field
            var id = $(select).attr("class").split(" ")[0], opval =  $(subSectionOption).attr("value");
            $("input#articleSubsection-" + id).attr("value",opval);

            //hide subsection error
            $("#articleSubsectionError-" + id).hide();

            if($(subSectionOption).val() > 0) {
                    jQuery.ajax({
                    url : "/author/manuscript/details/simpleSubsectionDescription",
                    dataType: "html",
                    data : {
                        "journalSubSectionId" : $(select).val()
                    },
                    success : function(data) {
                        $(desc).html(data);
                        $(desc).hide();
                        if($(desc).text().trim().length > 0) {
                            $(scopelink).show();
                        } else {
                            $(scopelink).hide();
                        }
                        },
                    error : function(request, errorType, errorThrown) {
                        $(desc).html('<h4 class="error">Sorry, an error occurred: ' + errorType + '</h4>');
                    },
                    cache : false,
                    timeout : 60000
                });
            } else {
                $(scopelink).hide();
            }
        });

       $('.showSectionDescription a').click(function(){
           $(this).siblings('.sectiondescription_result').toggle();
           $(this).toggleClass("viewScope");
       });

        $('.showSubSectionDescription a').click(function(){
            $(this).siblings('.subsectiondescription_result').toggle();
            $(this).toggleClass("viewScope");
        });

        //change subsection link (this will appear if the editor has suggested a subsection when reject and transfer)
        $('.changeSubsection').click(function(){
                    //hide subsection name
                    $(this).parent().children(".subSectionName").hide();

                    //hide view scope link
                    $(this).parent().children('.showSubSectionDescription').hide();

                    //hide the change link
                    $(this).hide();

                    //show section drop dawn
                    var subSectionSelect = $(this).parent().children('.journalSubSectionId');
                    subSectionSelect.show();

           //populating subsections for suggested (by the editor) section
           var sectionId = $(this).parents("tr").siblings(".section_tr").children(".section_td").children(".sectionName").attr("id");

            //setting hidden field value to 0
            var id = $(subSectionSelect).attr("class").split(" ")[0];
            $("input#articleSubsection-" + id).attr("value","0");

            jQuery.ajax({
                        url : "/author/manuscript/details/subSections",
                        dataType: "html",
                        data : {
                            "journalSectionId" : sectionId
                        },
                        success : function(data) {
                            $(subSectionSelect).html(data);
                        },
                        error : function(request, errorType, errorThrown) {
                            alert("Error while loading subsections")
                        },
                        cache : false,
                        timeout : 60000
            });
        });

        if($("span.sectionName").length > 0) {
            $("span.sectionName").each(function(){
                var manId =  $(this).attr("class").split(" ")[0];
                $("input#articleSection-" + manId).attr("value", $(this).attr("id"));
            });
        }

        if($("span.subSectionName").length > 0) {
            $("span.subSectionName").each(function(){
                var manId =  $(this).attr("class").split(" ")[0];
                $("input#articleSubsection-" + manId).attr("value", $(this).attr("id")); });
        }

        if($("select.journalSectionId").length > 0) {
            $("select.journalSectionId").each(function() {
                if($(this).is(":visible")) {
                    var o = $(this).children("option:selected");
                    var id = $(o).attr("class");
                    var optionValue =  $(o).attr("value");
                    $("input#articleSection-" + id).attr("value",optionValue);
                }
            });
        }

        //Article type
        $('.hideDetails').click(function(){

            // Show article type select
            $(this).parent().siblings('.article-type').show();

            // Hides man type name
            $(this).parent().siblings('.manuscriptTypeName').hide();

            // Hide that link and span
            $(this).parent().hide();
        });


        $(".article-type").change(function(){
            var o = $(this).children("option:selected"), id = $(o).attr("class"), opval =  $(o).attr("value");
            $("input#articleType-" + id).attr("value",opval);
            $("#articleTypeError-" + id).hide();
        });
         if($("span.manuscriptTypeName").length > 0) {
            $("span.manuscriptTypeName").each(function() { var cls = $(this).attr("class").split(" ")[0]; $("input#articleType-" + cls).attr("value", $(this).attr("id")); });
         }

        if($("select.article-type").length > 0) {
           $("select.article-type").each(function() { var o = $(this).children("option:selected"), id = $(o).attr("class"), opval =  $(o).attr("value"); $("input#articleType-" + id).attr("value",opval); });
        }
    }

    if(site.section == "my" && $("#inst-reports-accordion").length>0) {
        (function() {
            $("#inst-reports-accordion").click(function(e) {
                e.preventDefault();
                $("ul.secondary-nav li").removeClass("current");
                $(this).parents("li").toggleClass("current");
                $(this).siblings("ul.tertiary-nav").toggle();
            });
        })();
    };

    // FAQ / browse articles scroll
    jQuery('p.options a.options-toggler').each(function () {
        jQuery(this).removeClass('active').parents('.options').eq(0).next().hide()
    });


    if (jQuery('a[rel="photos_group"]').length != 0 || jQuery('a.lightbox').length != 0) {
        site.initializeLightbox();
    }
         if ($("a#whatismyip").length != 0) {
        site.initializeLightbox();
    }

    if ($("#flowplayer").offset()) {
        site.videoPlaylist(video_args);
    }

    if ($(".myPlayer").offset()) {
        site.singleVideo(".myPlayer");
    }

    if (site.page == "submission") {
     $("#subsection_dropdown").hide();

     if ($('#journalSubSectionId option').size() <= 1) {

             $("#journalSubSectionId").hide();
     } else {

          $("#journalSubSectionId").show();
     }

        site.loadJournalSection();
        site.loadJournalSubSection();
        site.loadJournalSubSectionDescription();
        site.loadSubjects();
    }
    if ((site.id == "3009" || site.id == "3008") && site.page == "home" || site.page == "collections" && window.bmcIsMobile != "mobile") {
        flowplayer("audio", "http://releases.flowplayer.org/swf/flowplayer-3.2.7.swf", {
            clip: {autoPlay: false},
            plugins: {
                controls: {fullscreen: false, autoHide: false, mute: false, volume: true}
            }
        });        
	if ((site.id == "3009") && site.page == "home" && window.bmcIsMobile != "mobile"){	
		flowplayer("audio2", "http://releases.flowplayer.org/swf/flowplayer-3.2.7.swf", {
            clip: {autoPlay: false},
            plugins: {
                controls: {fullscreen: false, autoHide: false, mute: false, volume: true}
            }
        });
    }
	}
    /* render video article flowplayer */
    if($("#videoarticle").length > 0) {
        $("#videoarticle").click(function(e) {
            e.preventDefault();
               /*$f("videoarticle", "http://bmcimg.com/flash/flowplayer/flowplayer-3.2.7.swf", {
                 buffering: true,
                 clip: {
                   autoPlay: true,
                   autoBuffering: true,
                   bufferLength: '3',
                   scaling: 'fit'
                 }
               }).ipad();*/
                $("#videoarticle").parent().html("<video controls='controls' width='450' height='285' preload='metadata' autoplay><source src='"+$(this).attr("href")+"' type='video/mp4'></source></video>");
               });
           }

    /* pop-up window for tracleer ad */
    if($('a[href*="Tracleer"]').length > 0) {
        $('a[href*="Tracleer"]').click(function(e) {
            e.preventDefault();
            window.open("http://"+site.url+"/advertising/tracleerpi.html","Tracleer","location=0,status=0,scrollbars=1, width=600,height=800", false)
        });
    }
    if(site.id=="9015" || site.id=="10324") {
        flowplayer("elifeplayer", "http://releases.flowplayer.org/swf/flowplayer-3.2.7.swf");
    }

    /* editor assignment radio-tr toggle */
     if(site.page == "submission" && $("table.moreinfo").length > 0)   {
     $("input[type=radio]").click(function() {
     $("table.moreinfo").find("tr.odd, tr.even").css("background-color","#eee");
     $(this).parents("tr").first().css("background-color","darkgrey");
    });
    }

    //Altmetric score
    if ($("#altmetric-badge-container").length > 0) {
        var getDOI = jQuery("#altmetric-badge-container").attr("data-doi");
        jQuery.ajax({
            type: "GET",
            url: "http://api.altmetric.com/v1/doi/"+getDOI,
            beforeSend: function(x) {
                if(x && x.overrideMimeType) {
                    x.overrideMimeType("application/j-son;charset=UTF-8");
                }
            },
            dataType: "jsonp",
            success: function(data){
                if(data.score>0) {
                    jQuery("#altmetric-badge-container").html("<a href='"+data.details_url+"'><img src='"+data.images.small+"'/></a><p>Altmetric score from <a href='http://altmetric.com/index.php'>Altmetric.com</a></p>");
                }
            }
        }); 
    }
    
    if (site.page == "china-microsite") {
        /*China gateways signup box  */
        if($("#signup-china").length > 0) {
            $("#signup-china").submit(function(e) {
                e.preventDefault();
                site.submitnewsletter({address: $(this).find("input#email").attr("value"), newsletters: "bmcupdates", action:true, messagecontainer: "#message-box"} );
            });
        }
        if($("#signup-china-news").length > 0) {
            $("#signup-china-news").submit(function(e) {
                e.preventDefault();
                site.submitnewsletter({address: $(this).find("input#email-news").attr("value"), newsletters: "bmcupdates", action:true, messagecontainer: "#message-box-news"} );
            });
        }
    }
	if (site.id == "9001" || site.id == "9013" || site.id == "9014") {
		$('#supplements_link').click(function(){
			$('#article_collection_and_supplements_panel .nav').find('.current').removeClass('current')
			$(this).parent().addClass('current');
			$('#supplements_content').show();
			$('#article_collection_content').hide();
			return false;
		});

		$('#article_collection_link').click(function(){
			$('#article_collection_and_supplements_panel .nav').find('.current').removeClass('current')
			$(this).parent().addClass('current');
			$('#article_collection_content').show();
			$('#supplements_content').hide();
			return false;
		});
	}

/*copied from mobile.js*/	
	site.allCookies = {
		getItem: function (sKey) {
			if (!sKey || !this.hasItem(sKey)) { return null; }
			return unescape(document.cookie.replace(new RegExp("(?:^|.*;\\s*)" + escape(sKey).replace(/[\-\.\+\*]/g, "\\$&") + "\\s*\\=\\s*((?:[^;](?!;))*[^;]?).*"), "$1"));
		},
		setItem: function (sKey, sValue, vEnd, sPath, sDomain, bSecure) {
			if (!sKey || /^(?:expires|max\-age|path|domain|secure)$/i.test(sKey)) { return; }
			var sExpires = "";
			if (vEnd) {
				switch (vEnd.constructor) {
					case Number:
						sExpires = vEnd === Infinity ? "; expires=Tue, 19 Jan 2038 03:14:07 GMT" : "; max-age=" + vEnd;
						break;
					case String:
						sExpires = "; expires=" + vEnd;
						break;
					case Date:
						sExpires = "; expires=" + vEnd.toGMTString();
						break;
				}
			}
			document.cookie = escape(sKey) + "=" + escape(sValue) + sExpires + (sDomain ? "; domain=" + sDomain : "") + (sPath ? "; path=" + sPath : "") + (bSecure ? "; secure" : "");
		},
	
		removeItem: function (sKey, sPath) {
		if (!sKey || !this.hasItem(sKey)) { return; }
		document.cookie = escape(sKey) + "=; expires=Thu, 01 Jan 1970 00:00:00 GMT" + (sPath ? "; path=" + sPath : "");
		},
	
		hasItem: function (sKey) {
		return (new RegExp("(?:^|;\\s*)" + escape(sKey).replace(/[\-\.\+\*]/g, "\\$&") + "\\s*\\=")).test(document.cookie);
		},
	
		keys: /* optional method: you can safely remove it! */ function () {
		var aKeys = document.cookie.replace(/((?:^|\s*;)[^\=]+)(?=;|$)|^\s*|\s*(?:\=[^;]*)?(?:\1|$)/g, "").split(/\s*(?:\=[^;]*)?;\s*/);
		for (var nIdx = 0; nIdx < aKeys.length; nIdx++) { aKeys[nIdx] = unescape(aKeys[nIdx]); }
		return aKeys;
		}
	};
			
	site.viewtype = function() {
        
            
            var queryString=location.search;
            var params=queryString.substring(1).split('&');
            for(var i=0; i<params.length; i++) {
                var pair=params[i].split('=');
                var key=decodeURIComponent(pair[0]);
                var dec=decodeURIComponent(pair[1]);
                if(key == "fmt_view" && dec == "classic") {
                    $("body").removeClass("mobile").addClass("classic");
                } 
                if(key == "fmt_view") {
                   // site.allCookies.setItem("fmt_view",dec,150000, "/");
                }
            }


        
            // use varnish method:
            // for varnish we are creatign the cookie before the page reloads so the mobile cookie is sent as part of the request header
           
            $("a#classic-view, a#mobile-view").click(function(e) {
                if($(this).attr("id") == "mobile-view") {
                    site.allCookies.setItem("fmt_view","mobile",150000, "/");
                } else if($(this).attr("id") == "classic-view") {
                     site.allCookies.setItem("fmt_view","classic",150000, "/");
                }
            });

            
        
    }
    
    //china gateway iframe hack for journal a-z
	
    if(site.page == "china-microsite" && jQuery("#chinaGatewayJournals").length > 0) {
        jQuery("#chinaGatewayJournals").load(function() {
            var chinaGatewayJournalsIframe = jQuery("#chinaGatewayJournals").contents();
            var browseBySubject = chinaGatewayJournalsIframe.find("#browse-by-subject-link");
            var journalAz = chinaGatewayJournalsIframe.find("#browse-by-az-link");
            if(endsWith(document.URL, '_en')) {
                chinaGatewayJournalsIframe.find('#journals li h2.toggler').addClass('active selected').next().show(); 
                browseBySubject.attr("href", browseBySubject.attr('href')+"_en");
                journalAz.attr("href", journalAz.attr('href')+"_en");
            }
            if(journalAz.parent().hasClass("current")) {
                chinaGatewayJournalsIframe.find("#journals-az a").click(function(e){
                    e.preventDefault();
                    var anchor = $(this).text().replace(/ /g,'');
                    $(document).scrollTop(chinaGatewayJournalsIframe.find("a[name='"+anchor+"']").offset().top +175);
                });
            }
            chinaGatewayJournalsIframe.find("li.align-right a").click(function(e){
                e.preventDefault();
                window.location.hash="content";
            });
            setIframeHeight(document.getElementById('chinaGatewayJournals'));
        });
    }
    
    if ($("#twitter_panel").length > 0) {
        // For journal 10172
        site.updateTwitterPanel();
    }
    if ($("a.flowplayer").length > 0) {
        // setup player
        site.setUpPlayer();
    }


site.renderTags = function(params) {

    var list = "",  array1 = [], array2 = [];

// creating two arrays to merge
    for (var i in params.array) {
        array1.push([i,parseInt(params.array[i])]);
    }

// merge the arrays and remove duplicates
    var array3 = arrayUnique(array1.concat(array2));

    // sort array and merge dupes
    function arrayUnique(array) {
        var a = array.concat();
        for(var i=0; i<a.length; ++i) {
            for(var j=i+1; j<a.length; ++j) {
                if(a[i][0].toLowerCase() == a[j][0].toLowerCase()) {
                    a.splice(j--, 1);
                    a[i][1] += a[j][1];
                }
            }
        }

        return a;
    };

    // sort alphabetically
    var sorted =array3.sort(function(a, b){
        if(a[0].toLowerCase() < b[0].toLowerCase()) return -1;
        if(a[0].toLowerCase() > b[0].toLowerCase()) return 1;
        return 0;
    });

    // print
    for(var i = 0; i < sorted.length; i++) {
        if(sorted[i]) {
            list = '<li value="' + sorted[i][1] + '"><a href="' + window.location.href  + "/" + sorted[i][0] + '">' + sorted[i][0].toLowerCase() + '</a></li>\n';
            $(params.container).append(list);
        }

    }

    // generate cloud
    $(params.container).tagcloud({type:"list",sizemin:11, sizemax:23, colormin:"#4760b0", colormax:"#4760b0", height: 300}).find("li");

}


site.trendingArticles = function(options) {

    window.alt_api = {results: [], articles: [], items: 0};

    // fetch altmetric results 
    $.ajax({url:options.uri, timeout: 30000, type: "GET", dataType: "jsonp", cache: true})
    .error(function(event, request, settings){ 

        $( "#error-message" ).append("<div class='msg msg-err'><div class='msg-err-img'></div><div class='msg-text'>There is currently no trending data available from Altmetric for this institution. Please try again later. </div>");
        $( "#loading" ).hide();
        // settings is error message
       // console.log(settings);

    })
    .success(function(data, textStatus, jqXHR) { 

        // console.log("succesfully fetched altmetric data");

        // if altmetric results are fewer then requested, set the counter accordingly
        if(data.results.length < options.size) { options.size = data.results.length; }

        if(data.results.length > 0) {

            var bmcsearch = [];

            for(var i = 0; i < options.size; i++) {

              // increment altmetric results
              window.alt_api.results.push(data.results[i]);

              bmcsearch[i] = {name : "doi", value : data.results[i].doi};

            }


             jso_configure({
                   "biomedcentral": {
                       client_id: "script-client",
                       redirect_uri: window.location,
                       authorization: "/oauth/authorize",
                       debug: true
                   }
               });

               $.oajax({
                   url: "/api/articles/find",
                   data: bmcsearch,
                   jso_provider: "biomedcentral",
                   jso_allowia: true,
                   dataType: "json",
                   cache: true,
                   type: "GET",
                   timeout: 30000,
                   success: function(data) {

                   // console.log("succesfully fetched bmc data");

                    for(var i = 0; i <= window.alt_api.results.length; i++) {
                        for(var j = 0; j <= data.length; j++) {
                            if(typeof window.alt_api.results[i] != "undefined" && typeof data[j] != "undefined" && window.alt_api.results[i].doi == data[j].doi) {
                                render(data[j]);
                            }
                        }
                    }

                    // cache the results
                    if (Modernizr.localstorage) {
                      var bmcaltresults = $("#trending_articles").html();
                      localStorage.setItem("bmcaltresults", bmcaltresults);
                      localStorage.setItem("inst", window.inst);
                    } 

                    $( "#loading" ).hide();

                   },
                   error: function(xhr, status) {
                      // console.log("There was an error: " + status);
                       $( "#error-message" ).append("<div class='msg msg-err'><div class='msg-err-img'></div><div class='msg-text'>Error requesting data from BioMed Central. Please try again later. </div>");
                       $( "#loading" ).hide();
                   }
               });

    } else {

        $("#no-trending").show();
        $( "#loading" ).hide();

    }
        
    });

    function render(entries) {

            function formatDate(datestring) {
                // format the date to bmc style
                var d = new XDate(datestring);
                return d.toString("d MMMM yyyy");

            }

            function renderAuthors(authorobj) {

                var authors = authorobj.slice();
                var tx = "";
               // console.log(authors);
                for(var i = 0; i < authors.length; i++) {


                if(authors.length == 1) {
                    tx+= authors[i].name;
                } else if(i == (authors.length-1)) {
                    tx +=  " and " + authors[i].name;
                } else if (i == (authors.length-2)) {
                    tx += authors[i].name;
                } else {
                    tx += authors[i].name + ", ";
                }
                

            }

            return tx;
        }
            


            var html = '';

                //for(var i = 0; i <= entries.length; i++) {

                     html = '<li class="articles-feed">';
                     html += '<div data-badge-type="donut" class="altmetric-embed" data-badge-popover="right" data-doi="' + entries.doi + '"></div>';
                     html += '<div class="article-citation"><p class="article_type">' + (entries["type"]) + '</p><blockquote><h4>' + '<a class="article_title" href="' + entries.url + '">';
                     html += entries.title + '</a></h4></blockquote>';
                     html += '<p class="citation">' + (renderAuthors(entries.authors));
                     html += '<span class="article-citation"><em>' + entries.journal["name"] + '</em> ' + entries.year.toString() + '<strong>, ' + entries.volume + '</strong>:' + entries.firstPage + ' (' + formatDate(entries.publicationDate) + ') </span></p>'; 
                     html += '<a href="' + entries.url + '/about">Accesses</a> | ' + '<a href="' + entries.url + '/about#citations">Citations</a>';
                     html += '</div></li>'; 

                     $(html).appendTo(options.container);



                    // run altmetric badge code
                    _altmetric_embed_init();

                //}

    }

    function errorHandler() {

        // no results so hide error and show message
        if($("#trending_articles li").length == 0) {

            $("#no-trending").show();
            $( "#error-message" ).hide();

        }

    }

}





$(document).ready(function() {
    if(site.page == "Trending") {

        if($("#institute").length > 0 && typeof window.inst == "string") {

            try {

                    window.inst = window.inst.replace(". ", ".");
               
                    site.trendingArticles({
                        uri: "http://api.altmetric.com/v1/citations/at?doi_prefix=10.1186&pubmed_q=" + encodeURIComponent(window.inst),
                        size: 10,
                        container: "#trending_articles"
                    });
                
            } catch(error) 
            {
               // console.log(error);
               //$( "#error-message" ).append("<div class='msg msg-err'><div class='msg-err-img'></div><div class='msg-text'>Error requesting data. Please try again later. </div>");

            }

        }

    }
	
	site.viewtype();
                 
 });
  
// Added for internal page link on Editorial Policies page for Policies FAQ tab to open- About editors section to open a sliding panel on same page.
 $('.FAQLink').click(function(){
        var ref = $('#Editors_pages_content li.tab_2').attr('class');
        $('.'+ref).find('h2.toggler').addClass('active').addClass('selected');
        $('.'+ref).find('ul.hidebeforeload').slideDown(300);

          $('#editors-navigation li.current li.tab_2').find('a').each(function(){
             $(this).removeClass('active');
             $(this).addClass('selectedLink');
     });

    });
});

// biome widget should be hidden if empty

function biomeHide() {
	if($("iframe#biome").length > 0) {
		if($("iframe#biome").attr("src").match("loaded") == null) {
			$("iframe#biome").detach();
		}
	}
}

$(document).ready(function() {
	var biomeTimer = window.setTimeout(function() {
	    biomeHide();
	}, 3000);
});


 
jQuery(window).load(function() {
    if($(".comments_feed1").length > 0) {
                 $(".quotes").delegate("li.comment a.comment-blurb-area", "click", function(e) {        
                    if($(this).hasClass("comment-blurb-gradient")) {
                        $(this).removeClass("comment-blurb-gradient").children().removeClass("down").addClass("up").parent().prev().css("height", "auto");
                    } else { 
                        $(this).addClass("comment-blurb-gradient").children().removeClass("up").addClass("down").parent().prev().css("height", "80px");   
                    }
                    e.preventDefault();
                });
            }

            $(".comments_feed1 li.comment").each(function() {

                if($(this).text().length < 300) {
                    $(this).div().css("height", "auto").next('a.comment-blurb-area').hide();
                }

 });
     
// Added for subnavigation in About editors section to open a sliding panel on same page.
   
$('#editors-navigation li.current ul.sub li').click(function(){
        var anchrRef = $(this).attr('class');
        var expandCollapseRef = $('#Editors_pages_content').find('li.'+anchrRef).attr('class');
       
        if(anchrRef == expandCollapseRef){
            resetToggle();
           
            $(this).removeClass('active');
            $(this).find('a').addClass('selectedLink').attr('disabled');
      
            $('#Editors_pages_content').find('li.'+anchrRef).find('h2.toggler').addClass('active').addClass('selected');
             $('#Editors_pages_content').find('li.'+anchrRef).find('ul.hidebeforeload:first').slideDown(300);
           // $(this).find('ul.hidebeforeload').slideDown(300);
        }
        else{
                resetToggle();
                $(this).removeClass('active');
                $(this).find('a').addClass('selectedLink').attr('disabled');
                                     
        }
            });



function resetToggle(){
     $('#Editors_pages_content li h2').removeClass('active').removeClass('selected');    
    
      $('#Editors_pages_content li ul.hidebeforeload').slideUp(200);
      $('#editors-navigation li.current li').find('a').each(function(){
             $(this).removeClass('selectedLink');
             $(this).addClass('active');
     });
}
    var loc = window.location.href;
    var matches = loc.match(/#(.*)#/);
    if(matches != null && matches.length == 2)
    {
        if(matches[1].indexOf('=') > -1)
        {
            var elem = matches[1].split('=')[1];       
            var anchorRef = $('.about-nav').find('a[href*='+elem+']');
            $(anchorRef).parent().trigger('click');     
        }
   }

   });
